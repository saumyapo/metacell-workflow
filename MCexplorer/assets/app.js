// MCexplorer: browser-side helpers

// selectize menus open on top of everything (not clipped by cards or panels)
(function setSelectizeDefaults() {
  if (window.Selectize && Selectize.defaults) { Selectize.defaults.dropdownParent = 'body'; return; }
  setTimeout(setSelectizeDefaults, 10);
})();
if (!window.EyeDropper) document.documentElement.classList.add('mcx-no-eyedropper');

(function wait() {
  if (!(window.Shiny && Shiny.addCustomMessageHandler && Shiny.InputBinding)) { setTimeout(wait, 30); return; }

  // disable inputs; the reason is shown on hover so it is clear how to turn them back on
  Shiny.addCustomMessageHandler('mcx-disable', function(m) {
    (m.ids || []).forEach(function(id) {
      var el = document.getElementById(id); if (!el) return;
      var grp = $(el).closest('.shiny-input-container, .form-group, .bslib-input-switch, .form-switch, .btn');
      if (!grp.length) grp = $(el);
      grp.toggleClass('mcx-disabled', !!m.disable);
      if (m.disable && m.reason) grp.attr('data-tip', m.reason); else grp.removeAttr('data-tip');
      $(el).find('input,select,button,textarea').addBack().prop('disabled', !!m.disable);
      if (el.selectize) { m.disable ? el.selectize.disable() : el.selectize.enable(); }
    });
  });
  Shiny.addCustomMessageHandler('mcx-disable-choice', function(m) {
    var inp = $('input[name="' + m.id + '"][value="' + m.value + '"]'), lab = inp.closest('label');
    inp.prop('disabled', !!m.disable); lab.toggleClass('mcx-off', !!m.disable);
    var base = lab.attr('data-base-tip');
    if (m.disable && m.reason) lab.attr('data-tip', m.reason); else if (base) lab.attr('data-tip', base); else lab.removeAttr('data-tip');
  });
  Shiny.addCustomMessageHandler('mcx-clear', function(id) { $('#' + id).val('').trigger('change'); });
  // show or hide one panel of a collapsible group, and open it
  Shiny.addCustomMessageHandler('mcx-acc', function(m) {
    var item = $('#' + m.acc + ' > .accordion-item[data-value="' + m.panel + '"]'); if (!item.length) return;
    item.toggleClass('mcx-hidden', !m.show);
    if (m.show && m.open) { var c = item.find('> .accordion-collapse'); if (!c.hasClass('show')) bootstrap.Collapse.getOrCreateInstance(c[0]).show(); }
  });
  Shiny.addCustomMessageHandler('mcx-datalist', function(m) {
    var dl = document.getElementById(m.id); if (!dl) return;
    dl.innerHTML = '';
    (m.values || []).forEach(function(v) { var o = document.createElement('option'); o.value = v; dl.appendChild(o); });
  });

  // draggable, tickable list: value = ticked items in display order
  function rankCount(el) {
    var all = el.querySelectorAll('.mcx-rank-item').length, on = $(el).find('.mcx-rank-item input:checked').length;
    el.querySelector('.mcx-rank-count').textContent = all ? (on + ' of ' + all + ' shown') : '';
  }
  var rank = new Shiny.InputBinding();
  $.extend(rank, {
    find: function(scope) { return $(scope).find('.mcx-rank'); },
    getValue: function(el) {
      return $(el).find('.mcx-rank-item').filter(function() { return $(this).find('input').prop('checked'); })
        .map(function() { return this.getAttribute('data-value'); }).get();
    },
    subscribe: function(el, cb) { $(el).on('mcxrank.mcx', function() { cb(false); }); },
    unsubscribe: function(el) { $(el).off('.mcx'); },
    getRatePolicy: function() { return { policy: 'debounce', delay: 300 }; },
    initialize: function(el) {
      var list = el.querySelector('.mcx-rank-list');
      Sortable.create(list, { animation: 120, ghostClass: 'mcx-rank-ghost', filter: 'input', preventOnFilter: false,
                              onEnd: function() { $(el).trigger('mcxrank'); } });
      $(el).on('change', '.mcx-rank-item input', function(e) {
        e.stopPropagation(); $(this).closest('.mcx-rank-item').toggleClass('off', !this.checked); rankCount(el); $(el).trigger('mcxrank');
      });
      $(el).on('click', '.mcx-rank-lbl, .mcx-rank-sw', function() {
        var cb = $(this).closest('.mcx-rank-item').find('input'); cb.prop('checked', !cb.prop('checked')).trigger('change');
      });
      $(el).on('click', '.mcx-rank-all, .mcx-rank-none', function(e) {
        e.preventDefault(); var on = $(this).hasClass('mcx-rank-all');
        $(el).find('.mcx-rank-item input').prop('checked', on); $(el).find('.mcx-rank-item').toggleClass('off', !on);
        rankCount(el); $(el).trigger('mcxrank');
      });
    },
    receiveMessage: function(el, data) {
      var list = el.querySelector('.mcx-rank-list'); list.innerHTML = '';
      (data.items || []).forEach(function(it) {
        var row = document.createElement('div');
        row.className = 'mcx-rank-item' + (it.on ? '' : ' off'); row.setAttribute('data-value', it.value);
        row.innerHTML = '<span class="mcx-rank-handle" title="Drag to reorder"></span><input type="checkbox"' + (it.on ? ' checked' : '') +
                        '><span class="mcx-rank-sw"></span><span class="mcx-rank-lbl"></span>';
        row.querySelector('.mcx-rank-lbl').textContent = it.label;
        var sw = row.querySelector('.mcx-rank-sw');
        if (it.color && typeof it.color === 'string') sw.style.background = it.color; else sw.style.display = 'none';
        list.appendChild(row);
      });
      rankCount(el); $(el).trigger('mcxrank');
    }
  });
  Shiny.inputBindings.register(rank, 'mcx.rank');
})();

// colour wheel (native picker) and eyedropper next to every colour box
$(document).on('change', 'input.mcx-native', function() {
  Shiny.setInputValue('mcx_colour', { id: this.getAttribute('data-for'), value: this.value, t: Date.now() }, { priority: 'event' });
});
$(document).on('click', '.mcx-dropper', function() {
  var id = this.getAttribute('data-for'); if (!window.EyeDropper) return;
  new EyeDropper().open().then(function(r) {
    Shiny.setInputValue('mcx_colour', { id: id, value: r.sRGBHex, t: Date.now() }, { priority: 'event' });
  }).catch(function() {});
});
$(document).on('change', 'input.shiny-colour-input', function() {
  var n = document.querySelector('input.mcx-native[data-for="' + this.id + '"]');
  if (n && /^#[0-9a-fA-F]{6}$/.test(this.value)) n.value = this.value.toLowerCase();
});

// Enter in the name box applies the name
$(document).on('keydown', '#ann_name', function(e) {
  if (e.key === 'Enter') { e.preventDefault(); Shiny.setInputValue('ann_enter', { name: this.value, t: Date.now() }, { priority: 'event' }); }
});

// copy a gene list
$(document).on('click', '.mcx-copy', function() {
  var b = this, txt = b.getAttribute('data-copy');
  var done = function() { var old = b.innerHTML; b.textContent = 'Copied'; setTimeout(function() { b.innerHTML = old; }, 1400); };
  if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(txt).then(done);
  else { var t = document.createElement('textarea'); t.value = txt; document.body.appendChild(t); t.select(); document.execCommand('copy'); t.remove(); done(); }
});

// "How to use" links do not jump the page
$(document).on('click', 'a.mcx-howto', function(e) { e.preventDefault(); });

// "How to use" windows live at the end of the page so nothing covers them
$(function() { $('.mcx-howto-modal').appendTo(document.body); });

// names list on the Annotate page
$(document).on('click', '.mcx-name-row [data-act]', function() {
  Shiny.setInputValue('ann_act', { act: this.getAttribute('data-act'), name: $(this).closest('.mcx-name-row').attr('data-name'), t: Date.now() }, { priority: 'event' });
});

// small hover hints ([data-tip]) float above the whole page, so panels never hide them
(function() {
  var box = null;
  function show(el) {
    var t = el.getAttribute('data-tip'); if (!t) return;
    if (!box) { box = document.createElement('div'); box.className = 'mcx-floattip'; document.body.appendChild(box); }
    box.textContent = t; box.style.display = 'block';
    var r = el.getBoundingClientRect(), w = box.offsetWidth, h = box.offsetHeight;
    var x = Math.min(Math.max(6, r.left + r.width / 2 - w / 2), window.innerWidth - w - 6);
    var y = r.bottom + 6; if (y + h > window.innerHeight - 4) y = r.top - h - 6;
    box.style.left = x + 'px'; box.style.top = y + 'px';
  }
  $(document).on('mouseenter', '[data-tip]', function() { show(this); });
  $(document).on('mouseleave click', '[data-tip]', function() { if (box) box.style.display = 'none'; });
})();

// - and + buttons next to number boxes
$(document).on('click', '.mcx-step', function() {
  var inp = document.getElementById(this.getAttribute('data-for')); if (!inp) return;
  var st = parseFloat(inp.step) || 1, v = parseFloat(inp.value); if (isNaN(v)) v = 0;
  v = v + st * parseFloat(this.getAttribute('data-dir'));
  if (inp.min !== '' && v < parseFloat(inp.min)) v = parseFloat(inp.min);
  if (inp.max !== '' && v > parseFloat(inp.max)) v = parseFloat(inp.max);
  var dec = (String(st).split('.')[1] || '').length;
  inp.value = v.toFixed(dec); $(inp).trigger('change');
});

// the app name in the header goes to the Home page (it does not reload anything)
$(document).on('click', '#mcx_brand', function(e) { e.preventDefault(); var a = document.querySelector('.navbar a[data-value="home"]'); if (a) a.click(); window.scrollTo({ top: 0 }); });

// "back to top" button once the page is scrolled down
$(function() {
  var b = $('<button type="button" class="mcx-totop" title="Back to the top"><svg viewBox="0 0 16 16" width="14" height="14"><path d="M8 3l5 6H9v4H7V9H3z" fill="currentColor"/></svg></button>').appendTo(document.body);
  b.on('click', function() { window.scrollTo({ top: 0, behavior: 'smooth' }); });
  $(window).on('scroll', function() { b.toggleClass('show', window.scrollY > 300); });
});
// tab names as hover text, for narrow windows where only the icons show
// tab hover text: name and what the tab is for, in the app's own hint style
$(function() { $('.navbar .nav-link').each(function() { var n = $(this).find('.mcx-navtxt'); if (!n.length) return;
  this.removeAttribute('title'); this.setAttribute('data-tip', n.text() + ': ' + (n.attr('data-desc') || '')); }); });

// Annotate heatmap: the server sends where columns and blocks sit in the image (mcx-hm-geom) and which columns
// are selected (mcx-hm-sel). Hovering shows the block under the mouse, clicking or dragging sends hm_pick.
(function() {
  var G = null, SEL = [], drag = null;
  function wrap() { return document.getElementById('hm_wrap'); }
  function img() { var w = wrap(); return w ? w.querySelector('img') : null; }
  function colAt(px) { return Math.round(G.x0 + px / G.W * (G.x1 - G.x0)); }           // css px -> column
  function pxOf(c) { return (c - G.x0) / (G.x1 - G.x0) * G.W; }                        // column centre -> css px
  function runOf(c) { if (!G || !G.runs) return null; var lo = 0, hi = G.runs.length - 1;
    while (lo <= hi) { var m = (lo + hi) >> 1, r = G.runs[m]; if (c < r.s) hi = m - 1; else if (c > r.e) lo = m + 1; else return r; } return null; }
  function box(cls, s, e) { var d = document.createElement('div'); d.className = cls;
    d.style.left = pxOf(s - 0.5) + 'px'; d.style.width = Math.max(2, pxOf(e + 0.5) - pxOf(s - 0.5)) + 'px';
    d.style.top = G.top + 'px'; d.style.height = (G.bottom - G.top) + 'px'; return d; }
  function drawSel() { var w = wrap(); if (!w || !G) return; var L = w.querySelector('.mcx-hm-layer'); if (!L) return;
    L.querySelectorAll('.mcx-hm-sel').forEach(function(x) { x.remove(); });
    // zoomed in: only a thin outline, so the data stays readable; nothing when everything shown is selected
    if (G.zoom && SEL.length === 1 && SEL[0].s === 1 && SEL[0].e === G.nc) return;
    SEL.forEach(function(r) { if (r.e <= G.nc) L.appendChild(box('mcx-hm-sel' + (G.zoom ? ' outline' : ''), r.s, r.e)); }); }
  function pos(ev) { var im = img(); if (!im || !G) return null; var r = im.getBoundingClientRect();
    var x = ev.clientX - r.left, y = ev.clientY - r.top; var c = colAt(x);
    return { x: x, y: y, c: c, inside: c >= 1 && c <= G.nc && y >= G.top - 2 && y <= G.bottom + 2 }; }
  function hideHover() { var w = wrap(); if (!w) return; w.querySelectorAll('.mcx-hm-hover').forEach(function(x) { x.remove(); });
    var t = w.querySelector('.mcx-hm-tip'); if (t) t.style.display = 'none'; }
  Shiny.addCustomMessageHandler('mcx-hm-geom', function(g) { G = g; drawSel(); });
  Shiny.addCustomMessageHandler('mcx-hm-sel', function(m) { SEL = m.runs || []; drawSel(); });
  $(document).on('shiny:value', function(e) { if (e.name === 'hm_plot') setTimeout(drawSel, 0); });
  $(document).on('mousemove', '#hm_wrap', function(ev) {
    var p = pos(ev), w = wrap(); if (!p) return;
    w.style.cursor = p.inside ? (drag ? 'col-resize' : 'pointer') : 'default';
    hideHover(); if (!p.inside && !drag) return;
    var L = w.querySelector('.mcx-hm-layer'), t = w.querySelector('.mcx-hm-tip');
    if (drag) { var a = Math.max(1, Math.min(drag.c, p.c)), b = Math.min(G.nc, Math.max(drag.c, p.c));
      var ra = runOf(a), rb = runOf(b); if (ra && G.runs.length < G.nc) a = ra.s; if (rb && G.runs.length < G.nc) b = rb.e;
      L.appendChild(box('mcx-hm-hover drag', a, b)); return; }
    var r = runOf(p.c); if (!r) return;
    L.appendChild(box('mcx-hm-hover', r.s, r.e));
    t.innerHTML = r.info || ''; t.style.display = 'block';
    var tw = t.offsetWidth, x = Math.min(Math.max(4, p.x + 14), w.clientWidth - tw - 4);
    t.style.left = x + 'px'; t.style.top = Math.max(0, p.y - t.offsetHeight - 10) + 'px';
  });
  $(document).on('mouseleave', '#hm_wrap', function() { hideHover(); });
  $(document).on('mousedown', '#hm_wrap', function(ev) {
    if (ev.button !== 0) return; var p = pos(ev); if (!p || !p.inside) return;
    ev.preventDefault(); drag = { c: p.c, x: p.x, add: ev.shiftKey || ev.ctrlKey || ev.metaKey };
  });
  $(document).on('mouseup', function(ev) {
    if (!drag) return; var d = drag; drag = null; var p = pos(ev); hideHover(); if (!p) return;
    var c1 = Math.abs(p.x - d.x) < 4 ? d.c : Math.max(1, Math.min(G.nc, p.c));
    Shiny.setInputValue('hm_pick', { c0: d.c, c1: c1, add: d.add, t: Date.now() }, { priority: 'event' });
  });
  // double-click a block: zoom into it (double-click again to zoom back out)
  $(document).on('dblclick', '#hm_wrap', function(ev) { var p = pos(ev); if (!p) return;
    Shiny.setInputValue('hm_dbl', { c: p.inside ? p.c : -1, t: Date.now() }, { priority: 'event' }); });
  // a click anywhere outside the heatmap columns clears the selection, except on controls and the side panels
  var KEEP = '.mcx-insp, .sidebar, .modal, .mcx-blockbar, .navbar, .tooltip, .card-header, .card-footer, .mcx-fold, .dataTables_wrapper, ' +
             '.selectize-dropdown, button, input, select, textarea, a, label, .mcx-totop';
  $(document).on('mousedown', function(ev) {
    if (!G || !SEL.length || !$('#hm_wrap').is(':visible') || ev.button !== 0) return;
    if ($(ev.target).closest(KEEP).length) return;
    if ($(ev.target).closest('#hm_wrap').length) { var p = pos(ev); if (p && p.inside) return; }
    Shiny.setInputValue('sel_esc', Date.now(), { priority: 'event' });
  });
  $(document).on('keydown', function(ev) {
    if (ev.key === 'Escape' && $('#hm_wrap').is(':visible') && !$('.modal.show').length) Shiny.setInputValue('sel_esc', Date.now(), { priority: 'event' });
  });
})();

// tooltips never stay behind after a click
$(document).on('click', function() { $('.tooltip').remove(); });
