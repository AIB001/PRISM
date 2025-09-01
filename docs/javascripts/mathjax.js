// MathJax configuration for PRISM documentation
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
    tags: 'ams',
    packages: {
      '[+]': ['mathtools', 'cases', 'braket', 'mhchem']
    }
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  },
  loader: {
    load: ['[tex]/mathtools', '[tex]/cases', '[tex]/braket', '[tex]/mhchem']
  },
  chtml: {
    scale: 1.0,
    displayAlign: 'left',
    displayIndent: '2em'
  },
  svg: {
    scale: 1.0,
    displayAlign: 'left',
    displayIndent: '2em'
  }
};

// Custom PRISM documentation scripts
document$.subscribe(function() {
  // Add copy functionality to code blocks
  const codeBlocks = document.querySelectorAll('pre > code');
  codeBlocks.forEach(block => {
    // Add line numbers for long code blocks
    const lines = block.textContent.split('\n');
    if (lines.length > 5) {
      block.classList.add('line-numbers');
    }
  });

  // Smooth scroll for anchor links
  document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function (e) {
      e.preventDefault();
      const target = document.querySelector(this.getAttribute('href'));
      if (target) {
        target.scrollIntoView({
          behavior: 'smooth',
          block: 'start'
        });
      }
    });
  });

  // Add interactive tooltips
  const addTooltip = (element, text) => {
    element.setAttribute('data-tooltip', text);
    element.classList.add('has-tooltip');
  };

  // Add tooltips to abbreviations
  document.querySelectorAll('abbr').forEach(abbr => {
    if (abbr.title) {
      addTooltip(abbr, abbr.title);
    }
  });

  // Enhance tables with sorting capability
  const tables = document.querySelectorAll('table');
  tables.forEach(table => {
    // Add class for styling
    table.classList.add('sortable');
    
    // Add sort indicators to headers
    const headers = table.querySelectorAll('th');
    headers.forEach((header, index) => {
      header.style.cursor = 'pointer';
      header.setAttribute('data-sort-index', index);
      
      // Add sort icon
      const sortIcon = document.createElement('span');
      sortIcon.className = 'sort-icon';
      sortIcon.innerHTML = ' ⇅';
      header.appendChild(sortIcon);
      
      // Add click handler
      header.addEventListener('click', () => sortTable(table, index));
    });
  });

  // Progressive image loading
  const images = document.querySelectorAll('img[data-src]');
  const imageObserver = new IntersectionObserver((entries, observer) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        const img = entry.target;
        img.src = img.dataset.src;
        img.removeAttribute('data-src');
        imageObserver.unobserve(img);
      }
    });
  });

  images.forEach(img => imageObserver.observe(img));

  // Add animation to feature cards
  const cards = document.querySelectorAll('.grid.cards > *');
  const cardObserver = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        entry.target.classList.add('animate-in');
      }
    });
  }, { threshold: 0.1 });

  cards.forEach(card => cardObserver.observe(card));

  // Enhanced search experience
  const searchInput = document.querySelector('.md-search__input');
  if (searchInput) {
    // Add search suggestions
    searchInput.addEventListener('input', function() {
      const query = this.value.toLowerCase();
      // Implement search suggestions logic here
    });
  }

  // Add keyboard shortcuts
  document.addEventListener('keydown', function(e) {
    // Ctrl/Cmd + K for search
    if ((e.ctrlKey || e.metaKey) && e.key === 'k') {
      e.preventDefault();
      const searchInput = document.querySelector('.md-search__input');
      if (searchInput) {
        searchInput.focus();
      }
    }
    
    // Escape to close search
    if (e.key === 'Escape') {
      const searchInput = document.querySelector('.md-search__input');
      if (searchInput && searchInput === document.activeElement) {
        searchInput.blur();
      }
    }
  });

  // Add reading progress indicator
  const article = document.querySelector('article');
  if (article) {
    const progressBar = document.createElement('div');
    progressBar.className = 'reading-progress';
    progressBar.style.cssText = `
      position: fixed;
      top: 0;
      left: 0;
      height: 3px;
      background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
      z-index: 1000;
      transition: width 0.2s ease;
    `;
    document.body.appendChild(progressBar);

    window.addEventListener('scroll', () => {
      const scrollTop = window.scrollY;
      const scrollHeight = article.scrollHeight - window.innerHeight;
      const progress = (scrollTop / scrollHeight) * 100;
      progressBar.style.width = `${progress}%`;
    });
  }

  // Add code syntax highlighting enhancements
  if (typeof Prism !== 'undefined') {
    Prism.highlightAll();
  }

  // Initialize mermaid diagrams if present
  if (typeof mermaid !== 'undefined') {
    mermaid.initialize({
      startOnLoad: true,
      theme: 'default',
      themeVariables: {
        primaryColor: '#4051b5',
        primaryTextColor: '#fff',
        primaryBorderColor: '#7c4dff',
        lineColor: '#667eea',
        secondaryColor: '#764ba2',
        tertiaryColor: '#9575cd'
      }
    });
  }
});

// Table sorting function
function sortTable(table, columnIndex) {
  const tbody = table.querySelector('tbody');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  const header = table.querySelectorAll('th')[columnIndex];
  const isAscending = header.classList.contains('sort-asc');
  
  // Sort rows
  rows.sort((a, b) => {
    const aValue = a.cells[columnIndex].textContent.trim();
    const bValue = b.cells[columnIndex].textContent.trim();
    
    // Try to parse as numbers
    const aNum = parseFloat(aValue);
    const bNum = parseFloat(bValue);
    
    if (!isNaN(aNum) && !isNaN(bNum)) {
      return isAscending ? bNum - aNum : aNum - bNum;
    }
    
    // Sort as strings
    return isAscending 
      ? bValue.localeCompare(aValue)
      : aValue.localeCompare(bValue);
  });
  
  // Update table
  tbody.innerHTML = '';
  rows.forEach(row => tbody.appendChild(row));
  
  // Update sort indicators
  table.querySelectorAll('th').forEach(th => {
    th.classList.remove('sort-asc', 'sort-desc');
  });
  
  header.classList.add(isAscending ? 'sort-desc' : 'sort-asc');
}

// Add custom animations
const style = document.createElement('style');
style.textContent = `
  @keyframes fadeInUp {
    from {
      opacity: 0;
      transform: translateY(20px);
    }
    to {
      opacity: 1;
      transform: translateY(0);
    }
  }
  
  .animate-in {
    animation: fadeInUp 0.5s ease forwards;
  }
  
  .has-tooltip {
    position: relative;
    cursor: help;
  }
  
  .has-tooltip:hover::after {
    content: attr(data-tooltip);
    position: absolute;
    bottom: 100%;
    left: 50%;
    transform: translateX(-50%);
    padding: 0.5rem;
    background: rgba(0, 0, 0, 0.8);
    color: white;
    border-radius: 4px;
    font-size: 0.875rem;
    white-space: nowrap;
    z-index: 1000;
    pointer-events: none;
  }
  
  .sort-asc .sort-icon::after {
    content: ' ↑';
  }
  
  .sort-desc .sort-icon::after {
    content: ' ↓';
  }
`;
document.head.appendChild(style);