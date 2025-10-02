#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Real-time Monitoring Dashboard

Provides real-time monitoring, alerting, and dashboard capabilities
for PRISM molecular dynamics simulations with web interface and
automated notifications.
"""

import time
import json
import threading
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, Any, List, Optional, Callable, Union
from dataclasses import dataclass, asdict
from enum import Enum
import asyncio
import websockets
import logging


class AlertLevel(Enum):
    """Alert severity levels"""
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


@dataclass
class Alert:
    """Alert message"""
    timestamp: datetime
    level: AlertLevel
    title: str
    message: str
    source: str
    context: Dict[str, Any]
    resolved: bool = False
    
    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data['timestamp'] = self.timestamp.isoformat()
        data['level'] = self.level.value
        return data


@dataclass
class MonitoringRule:
    """Monitoring rule for automated alerts"""
    name: str
    description: str
    condition: Callable[[Dict[str, Any]], bool]
    alert_level: AlertLevel
    cooldown_minutes: int = 5
    last_triggered: Optional[datetime] = None
    
    def should_trigger(self, data: Dict[str, Any]) -> bool:
        """Check if rule should trigger"""
        if self.last_triggered:
            elapsed = datetime.now() - self.last_triggered
            if elapsed.total_seconds() < self.cooldown_minutes * 60:
                return False
        
        return self.condition(data)
    
    def trigger(self) -> None:
        """Mark rule as triggered"""
        self.last_triggered = datetime.now()


class SimulationMonitor:
    """Real-time simulation monitoring system"""
    
    def __init__(self, simulation_id: str):
        self.simulation_id = simulation_id
        self.alerts: List[Alert] = []
        self.monitoring_data: Dict[str, Any] = {}
        self.callbacks: List[Callable[[Dict[str, Any]], None]] = []
        self.rules: List[MonitoringRule] = []
        self.monitoring_active = False
        self.monitoring_thread = None
        self.websocket_server = None
        self.connected_clients = set()
        
        # Setup default monitoring rules
        self._setup_default_rules()
    
    def _setup_default_rules(self):
        """Setup default monitoring rules"""
        # High CPU usage rule
        self.add_rule(
            "high_cpu_usage",
            "Alert when CPU usage exceeds 95%",
            lambda data: data.get('cpu_percent', 0) > 95,
            AlertLevel.WARNING
        )
        
        # High memory usage rule
        self.add_rule(
            "high_memory_usage",
            "Alert when memory usage exceeds 16GB",
            lambda data: data.get('memory_mb', 0) > 16000,
            AlertLevel.WARNING
        )
        
        # Simulation failure rule
        self.add_rule(
            "simulation_failure",
            "Alert when simulation encounters errors",
            lambda data: data.get('errors', 0) > 0,
            AlertLevel.ERROR
        )
        
        # Poor performance rule
        self.add_rule(
            "poor_performance",
            "Alert when simulation performance is poor",
            lambda data: data.get('ns_per_day', float('inf')) < 0.5,
            AlertLevel.WARNING,
            cooldown_minutes=30
        )
        
        # Disk space rule
        self.add_rule(
            "low_disk_space",
            "Alert when available disk space is low",
            lambda data: data.get('disk_free_gb', float('inf')) < 10,
            AlertLevel.CRITICAL
        )
    
    def add_rule(self, name: str, description: str, condition: Callable,
                 alert_level: AlertLevel, cooldown_minutes: int = 5):
        """Add monitoring rule"""
        rule = MonitoringRule(
            name=name,
            description=description,
            condition=condition,
            alert_level=alert_level,
            cooldown_minutes=cooldown_minutes
        )
        self.rules.append(rule)
    
    def add_callback(self, callback: Callable[[Dict[str, Any]], None]):
        """Add callback for monitoring updates"""
        self.callbacks.append(callback)
    
    def update_monitoring_data(self, data: Dict[str, Any]):
        """Update monitoring data and check rules"""
        self.monitoring_data.update(data)
        self.monitoring_data['last_update'] = datetime.now()
        
        # Check monitoring rules
        for rule in self.rules:
            if rule.should_trigger(self.monitoring_data):
                alert = Alert(
                    timestamp=datetime.now(),
                    level=rule.alert_level,
                    title=rule.name.replace('_', ' ').title(),
                    message=rule.description,
                    source=f"monitoring_rule:{rule.name}",
                    context=self.monitoring_data.copy()
                )
                self.add_alert(alert)
                rule.trigger()
        
        # Notify callbacks
        for callback in self.callbacks:
            try:
                callback(self.monitoring_data)
            except Exception as e:
                logging.error(f"Callback error: {e}")
        
        # Broadcast to websocket clients
        asyncio.create_task(self._broadcast_update())
    
    def add_alert(self, alert: Alert):
        """Add new alert"""
        self.alerts.append(alert)
        
        # Keep only recent alerts (last 24 hours)
        cutoff = datetime.now() - timedelta(hours=24)
        self.alerts = [a for a in self.alerts if a.timestamp > cutoff]
        
        logging.info(f"Alert: {alert.title} - {alert.message}")
    
    def resolve_alert(self, alert_index: int):
        """Resolve an alert"""
        if 0 <= alert_index < len(self.alerts):
            self.alerts[alert_index].resolved = True
    
    def get_dashboard_data(self) -> Dict[str, Any]:
        """Get current dashboard data"""
        current_alerts = [a for a in self.alerts if not a.resolved]
        
        return {
            'simulation_id': self.simulation_id,
            'timestamp': datetime.now().isoformat(),
            'monitoring_data': self.monitoring_data,
            'alerts': {
                'total': len(current_alerts),
                'critical': len([a for a in current_alerts if a.level == AlertLevel.CRITICAL]),
                'error': len([a for a in current_alerts if a.level == AlertLevel.ERROR]),
                'warning': len([a for a in current_alerts if a.level == AlertLevel.WARNING]),
                'recent': [a.to_dict() for a in current_alerts[-10:]]
            },
            'rules': [
                {
                    'name': rule.name,
                    'description': rule.description,
                    'level': rule.alert_level.value,
                    'last_triggered': rule.last_triggered.isoformat() if rule.last_triggered else None
                }
                for rule in self.rules
            ]
        }
    
    async def _broadcast_update(self):
        """Broadcast update to websocket clients"""
        if self.connected_clients:
            dashboard_data = self.get_dashboard_data()
            message = json.dumps(dashboard_data)
            
            # Send to all connected clients
            disconnected = set()
            for client in self.connected_clients:
                try:
                    await client.send(message)
                except websockets.exceptions.ConnectionClosed:
                    disconnected.add(client)
            
            # Remove disconnected clients
            self.connected_clients -= disconnected
    
    async def websocket_handler(self, websocket, path):
        """Handle websocket connections"""
        self.connected_clients.add(websocket)
        
        # Send initial data
        dashboard_data = self.get_dashboard_data()
        await websocket.send(json.dumps(dashboard_data))
        
        try:
            async for message in websocket:
                # Handle client messages if needed
                pass
        except websockets.exceptions.ConnectionClosed:
            pass
        finally:
            self.connected_clients.discard(websocket)
    
    def start_websocket_server(self, host: str = "localhost", port: int = 8765):
        """Start websocket server for real-time updates"""
        async def server():
            self.websocket_server = await websockets.serve(
                self.websocket_handler, host, port
            )
            await self.websocket_server.wait_closed()
        
        # Run in separate thread
        def run_server():
            asyncio.run(server())
        
        server_thread = threading.Thread(target=run_server, daemon=True)
        server_thread.start()
        
        logging.info(f"Websocket server started on ws://{host}:{port}")
    
    def stop_websocket_server(self):
        """Stop websocket server"""
        if self.websocket_server:
            self.websocket_server.close()


class HTMLDashboard:
    """HTML dashboard generator"""
    
    @staticmethod
    def generate_dashboard_html(monitor: SimulationMonitor, 
                              websocket_url: str = "ws://localhost:8765") -> str:
        """Generate HTML dashboard"""
        dashboard_data = monitor.get_dashboard_data()
        
        html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>PRISM Simulation Dashboard</title>
            <style>
                body {{
                    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                    margin: 0;
                    padding: 20px;
                    background-color: #f5f5f5;
                }}
                .header {{
                    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                    color: white;
                    padding: 20px;
                    border-radius: 10px;
                    margin-bottom: 20px;
                }}
                .dashboard {{
                    display: grid;
                    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
                    gap: 20px;
                }}
                .card {{
                    background: white;
                    border-radius: 10px;
                    padding: 20px;
                    box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
                }}
                .metric {{
                    display: flex;
                    justify-content: space-between;
                    padding: 10px 0;
                    border-bottom: 1px solid #eee;
                }}
                .metric:last-child {{
                    border-bottom: none;
                }}
                .metric-value {{
                    font-weight: bold;
                    color: #333;
                }}
                .alert {{
                    padding: 10px;
                    margin: 5px 0;
                    border-radius: 5px;
                    border-left: 4px solid;
                }}
                .alert.critical {{
                    background-color: #ffeaea;
                    border-color: #f44336;
                }}
                .alert.error {{
                    background-color: #fff3e0;
                    border-color: #ff9800;
                }}
                .alert.warning {{
                    background-color: #fff8e1;
                    border-color: #ffc107;
                }}
                .alert.info {{
                    background-color: #e3f2fd;
                    border-color: #2196f3;
                }}
                .status {{
                    padding: 5px 10px;
                    border-radius: 15px;
                    font-size: 12px;
                    font-weight: bold;
                    text-transform: uppercase;
                }}
                .status.connected {{
                    background-color: #4caf50;
                    color: white;
                }}
                .status.disconnected {{
                    background-color: #f44336;
                    color: white;
                }}
                .chart-container {{
                    height: 200px;
                    margin-top: 20px;
                }}
                #connectionStatus {{
                    position: fixed;
                    top: 20px;
                    right: 20px;
                    z-index: 1000;
                }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>PRISM Simulation Dashboard</h1>
                <p>Simulation ID: {dashboard_data['simulation_id']}</p>
                <div id="connectionStatus">
                    <span class="status disconnected">Connecting...</span>
                </div>
            </div>
            
            <div class="dashboard">
                <div class="card">
                    <h3>System Metrics</h3>
                    <div id="systemMetrics">
                        <div class="metric">
                            <span>CPU Usage</span>
                            <span class="metric-value" id="cpuUsage">-</span>
                        </div>
                        <div class="metric">
                            <span>Memory Usage</span>
                            <span class="metric-value" id="memoryUsage">-</span>
                        </div>
                        <div class="metric">
                            <span>Performance</span>
                            <span class="metric-value" id="performance">-</span>
                        </div>
                        <div class="metric">
                            <span>Last Update</span>
                            <span class="metric-value" id="lastUpdate">-</span>
                        </div>
                    </div>
                </div>
                
                <div class="card">
                    <h3>Alert Summary</h3>
                    <div id="alertSummary">
                        <div class="metric">
                            <span>Critical</span>
                            <span class="metric-value" id="criticalAlerts">0</span>
                        </div>
                        <div class="metric">
                            <span>Errors</span>
                            <span class="metric-value" id="errorAlerts">0</span>
                        </div>
                        <div class="metric">
                            <span>Warnings</span>
                            <span class="metric-value" id="warningAlerts">0</span>
                        </div>
                        <div class="metric">
                            <span>Total Active</span>
                            <span class="metric-value" id="totalAlerts">0</span>
                        </div>
                    </div>
                </div>
                
                <div class="card">
                    <h3>Recent Alerts</h3>
                    <div id="recentAlerts">
                        <p>No alerts</p>
                    </div>
                </div>
                
                <div class="card">
                    <h3>Performance Chart</h3>
                    <div class="chart-container">
                        <canvas id="performanceChart" width="400" height="200"></canvas>
                    </div>
                </div>
            </div>
            
            <script>
                let socket;
                let performanceData = [];
                
                function connectWebSocket() {{
                    socket = new WebSocket('{websocket_url}');
                    
                    socket.onopen = function(event) {{
                        document.getElementById('connectionStatus').innerHTML = 
                            '<span class="status connected">Connected</span>';
                    }};
                    
                    socket.onmessage = function(event) {{
                        const data = JSON.parse(event.data);
                        updateDashboard(data);
                    }};
                    
                    socket.onclose = function(event) {{
                        document.getElementById('connectionStatus').innerHTML = 
                            '<span class="status disconnected">Disconnected</span>';
                        
                        // Retry connection after 5 seconds
                        setTimeout(connectWebSocket, 5000);
                    }};
                    
                    socket.onerror = function(error) {{
                        console.error('WebSocket error:', error);
                    }};
                }}
                
                function updateDashboard(data) {{
                    // Update system metrics
                    const monitoring = data.monitoring_data || {{}};
                    document.getElementById('cpuUsage').textContent = 
                        monitoring.cpu_percent ? monitoring.cpu_percent.toFixed(1) + '%' : '-';
                    document.getElementById('memoryUsage').textContent = 
                        monitoring.memory_mb ? (monitoring.memory_mb / 1024).toFixed(1) + ' GB' : '-';
                    document.getElementById('performance').textContent = 
                        monitoring.ns_per_day ? monitoring.ns_per_day.toFixed(1) + ' ns/day' : '-';
                    
                    if (monitoring.last_update) {{
                        const lastUpdate = new Date(monitoring.last_update);
                        document.getElementById('lastUpdate').textContent = lastUpdate.toLocaleTimeString();
                    }}
                    
                    // Update alert summary
                    const alerts = data.alerts || {{}};
                    document.getElementById('criticalAlerts').textContent = alerts.critical || 0;
                    document.getElementById('errorAlerts').textContent = alerts.error || 0;
                    document.getElementById('warningAlerts').textContent = alerts.warning || 0;
                    document.getElementById('totalAlerts').textContent = alerts.total || 0;
                    
                    // Update recent alerts
                    const recentAlertsDiv = document.getElementById('recentAlerts');
                    if (alerts.recent && alerts.recent.length > 0) {{
                        recentAlertsDiv.innerHTML = alerts.recent.map(alert => 
                            `<div class="alert ${{alert.level}}">
                                <strong>${{alert.title}}</strong><br>
                                ${{alert.message}}<br>
                                <small>${{new Date(alert.timestamp).toLocaleString()}}</small>
                            </div>`
                        ).join('');
                    }} else {{
                        recentAlertsDiv.innerHTML = '<p>No recent alerts</p>';
                    }}
                    
                    // Update performance chart
                    if (monitoring.ns_per_day !== undefined) {{
                        performanceData.push({{
                            timestamp: new Date(),
                            value: monitoring.ns_per_day
                        }});
                        
                        // Keep only last 50 data points
                        if (performanceData.length > 50) {{
                            performanceData.shift();
                        }}
                        
                        drawPerformanceChart();
                    }}
                }}
                
                function drawPerformanceChart() {{
                    const canvas = document.getElementById('performanceChart');
                    const ctx = canvas.getContext('2d');
                    
                    // Clear canvas
                    ctx.clearRect(0, 0, canvas.width, canvas.height);
                    
                    if (performanceData.length < 2) return;
                    
                    // Find min/max values
                    const values = performanceData.map(d => d.value);
                    const minValue = Math.min(...values);
                    const maxValue = Math.max(...values);
                    const range = maxValue - minValue || 1;
                    
                    // Draw axes
                    ctx.strokeStyle = '#ddd';
                    ctx.beginPath();
                    ctx.moveTo(40, 10);
                    ctx.lineTo(40, 190);
                    ctx.lineTo(canvas.width - 10, 190);
                    ctx.stroke();
                    
                    // Draw performance line
                    ctx.strokeStyle = '#2196f3';
                    ctx.lineWidth = 2;
                    ctx.beginPath();
                    
                    performanceData.forEach((point, index) => {{
                        const x = 40 + (index / (performanceData.length - 1)) * (canvas.width - 50);
                        const y = 190 - ((point.value - minValue) / range) * 180;
                        
                        if (index === 0) {{
                            ctx.moveTo(x, y);
                        }} else {{
                            ctx.lineTo(x, y);
                        }}
                    }});
                    
                    ctx.stroke();
                    
                    // Draw labels
                    ctx.fillStyle = '#666';
                    ctx.font = '12px Arial';
                    ctx.fillText('Performance (ns/day)', 10, 100);
                    ctx.fillText(maxValue.toFixed(1), 5, 15);
                    ctx.fillText(minValue.toFixed(1), 5, 195);
                }}
                
                // Connect to WebSocket on page load
                connectWebSocket();
            </script>
        </body>
        </html>
        """
        
        return html
    
    @staticmethod
    def save_dashboard(monitor: SimulationMonitor, output_file: Path,
                      websocket_url: str = "ws://localhost:8765"):
        """Save dashboard HTML to file"""
        html = HTMLDashboard.generate_dashboard_html(monitor, websocket_url)
        with open(output_file, 'w') as f:
            f.write(html)


# Convenience functions
def create_simulation_monitor(simulation_id: str) -> SimulationMonitor:
    """Create simulation monitor"""
    return SimulationMonitor(simulation_id)


def start_monitoring_dashboard(simulation_id: str, websocket_port: int = 8765,
                             save_html: bool = True) -> SimulationMonitor:
    """Start monitoring dashboard"""
    monitor = SimulationMonitor(simulation_id)
    
    # Start websocket server
    monitor.start_websocket_server(port=websocket_port)
    
    # Save HTML dashboard
    if save_html:
        dashboard_file = Path(f"dashboard_{simulation_id}.html")
        HTMLDashboard.save_dashboard(monitor, dashboard_file, f"ws://localhost:{websocket_port}")
        print(f"Dashboard saved to {dashboard_file}")
        print(f"Open {dashboard_file} in your browser to view the dashboard")
    
    return monitor