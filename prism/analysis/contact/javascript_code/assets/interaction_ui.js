/* ---------------------------------------------------------------------------
 * Interaction-report UI logic: theme toggle, occupancy slider, interaction-type
 * filters, legends and the interaction details table. All driven by the
 * INTERACTION_STYLE / INTERACTION_ORDER / RESIDUE_CLASS_COLORS / CONTACTS
 * globals injected by html_builder.py. No external dependencies.
 * ------------------------------------------------------------------------- */

function _dashToBorderStyle(dash) {
    if (!dash || dash.length === 0) return 'solid';
    // small repeating pattern reads as dotted, longer as dashed
    return (dash[0] <= 2) ? 'dotted' : 'dashed';
}

function _presentInteractionTypes() {
    // Types that actually occur in the data, ordered by INTERACTION_ORDER.
    const present = new Set();
    (typeof CONTACTS !== 'undefined' ? CONTACTS : []).forEach(c => {
        if (c.dominantType) present.add(c.dominantType);
        (c.interactionTypes || []).forEach(it => present.add(it.type));
    });
    const order = (typeof INTERACTION_ORDER !== 'undefined') ? INTERACTION_ORDER : [];
    const ordered = order.filter(t => present.has(t));
    // include any present-but-unordered types at the end
    present.forEach(t => { if (!ordered.includes(t)) ordered.push(t); });
    return ordered;
}

function toggleTheme() {
    const root = document.documentElement;
    const btn = document.getElementById('themeToggle');
    const dark = root.getAttribute('data-theme') === 'dark';
    root.setAttribute('data-theme', dark ? 'light' : 'dark');
    if (btn) btn.textContent = dark ? '🌙 Dark' : '☀️ Light';
}

function onOccupancyChange(value) {
    const v = parseFloat(value) / 100.0;
    const lbl = document.getElementById('occValue');
    if (lbl) lbl.textContent = Math.round(value) + '%';
    if (typeof contactMap !== 'undefined' && contactMap) {
        contactMap.minOccupancy = v;
    }
}

function _buildInteractionLegend() {
    const wrap = document.getElementById('interactionLegend');
    const empty = document.getElementById('interactionLegendEmpty');
    if (!wrap) return;
    const present = _presentInteractionTypes();
    if (present.length === 0) {
        wrap.innerHTML = '';
        if (empty) empty.style.display = 'block';
        return;
    }
    if (empty) empty.style.display = 'none';
    const summary = (typeof INTERACTION_SUMMARY !== 'undefined') ? INTERACTION_SUMMARY : {};
    wrap.innerHTML = present.map(t => {
        const s = INTERACTION_STYLE[t] || {label: t, color: '#999', dash: []};
        const style = _dashToBorderStyle(s.dash);
        const count = summary[t] !== undefined ? `<span class="count">${summary[t]} res</span>` : '';
        return `<div class="legend-item">
            <span class="legend-line" style="border-top-color:${s.color};border-top-style:${style};"></span>
            <span>${s.label}</span>${count}
        </div>`;
    }).join('');
}

function _buildResidueLegend() {
    const wrap = document.getElementById('residueLegend');
    if (!wrap || typeof RESIDUE_CLASS_COLORS === 'undefined') return;
    // only classes that appear among the contacts
    const present = new Set();
    (typeof CONTACTS !== 'undefined' ? CONTACTS : []).forEach(c => {
        if (c.residueClass) present.add(c.residueClass);
    });
    const labels = {
        hydrophobic: 'Hydrophobic', aromatic: 'Aromatic', positive: 'Positive (+)',
        negative: 'Negative (−)', polar: 'Polar', glycine: 'Glycine', other: 'Other'
    };
    const order = ['hydrophobic', 'aromatic', 'positive', 'negative', 'polar', 'glycine', 'other'];
    const shown = order.filter(c => present.has(c));
    wrap.innerHTML = (shown.length ? shown : order).map(c => {
        const col = RESIDUE_CLASS_COLORS[c] || '#888';
        return `<div class="legend-item">
            <span class="legend-dot" style="background:${col};"></span>
            <span>${labels[c] || c}</span>
        </div>`;
    }).join('');
}

function _buildTypeFilters() {
    const wrap = document.getElementById('typeFilters');
    if (!wrap) return;
    const present = _presentInteractionTypes();
    if (present.length === 0) { wrap.innerHTML = ''; return; }
    // start with all types active
    if (typeof contactMap !== 'undefined' && contactMap) {
        contactMap.activeTypes = new Set(present);
    }
    wrap.innerHTML = present.map(t => {
        const s = INTERACTION_STYLE[t] || {label: t, color: '#999', dash: []};
        const style = _dashToBorderStyle(s.dash);
        return `<span class="type-chip" data-type="${t}" onclick="toggleType(this)">
            <span class="swatch" style="border-top-color:${s.color};border-top-style:${style};"></span>${s.label}
        </span>`;
    }).join('');
}

function toggleType(el) {
    const t = el.getAttribute('data-type');
    if (typeof contactMap === 'undefined' || !contactMap) return;
    if (!contactMap.activeTypes) contactMap.activeTypes = new Set(_presentInteractionTypes());
    if (contactMap.activeTypes.has(t)) {
        contactMap.activeTypes.delete(t);
        el.classList.add('off');
    } else {
        contactMap.activeTypes.add(t);
        el.classList.remove('off');
    }
}

function _buildInteractionTable() {
    const body = document.getElementById('interactionTableBody');
    if (!body) return;
    const contacts = (typeof CONTACTS !== 'undefined' ? CONTACTS.slice() : []);
    const strength = c => (c.interactionTypes && c.interactionTypes.length)
        ? Math.max(c.frequency || 0, c.interactionTypes[0].occupancy)
        : (c.frequency || 0);
    contacts.sort((a, b) => strength(b) - strength(a));
    body.innerHTML = contacts.map(c => {
        const cls = c.residueClass || 'other';
        const clsCol = (typeof RESIDUE_CLASS_COLORS !== 'undefined' && RESIDUE_CLASS_COLORS[cls]) ? RESIDUE_CLASS_COLORS[cls] : '#888';
        let inter = '<span class="muted">contact only</span>';
        if (c.interactionTypes && c.interactionTypes.length) {
            inter = c.interactionTypes.map(it => {
                const s = (typeof INTERACTION_STYLE !== 'undefined') ? INTERACTION_STYLE[it.type] : null;
                const col = s ? s.color : '#888';
                const lbl = s ? s.label : it.type;
                return `<span class="pill" style="background:${col};" title="${(it.occupancy*100).toFixed(0)}% occupancy">${lbl} ${(it.occupancy*100).toFixed(0)}%</span>`;
            }).join(' ');
        }
        const top = (strength(c) * 100).toFixed(0);
        return `<tr>
            <td><b>${c.residue}</b></td>
            <td><span class="pill" style="background:${clsCol};">${cls}</span></td>
            <td>${inter}</td>
            <td><div class="occ-cell"><div class="occbar" style="width:${Math.max(6, top*0.6)}px;"></div><span>${top}%</span></div></td>
        </tr>`;
    }).join('');
}

window.addEventListener('DOMContentLoaded', () => {
    // Runs after utility_functions.js has created `contactMap`.
    try {
        _buildInteractionLegend();
        _buildResidueLegend();
        _buildTypeFilters();
        _buildInteractionTable();
    } catch (e) {
        console.error('Interaction UI init failed:', e);
    }
});
