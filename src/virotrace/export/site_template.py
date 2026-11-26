"""The explorer page: static HTML, CSS and JavaScript, no external requests.

Kept as a Python string so it ships with the package without package-data configuration.
``__TITLE__`` and ``__DATA__`` are replaced by ``site.render``.
"""

PAGE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root{
  --navy:#101D2D; --navy-2:#173B4D; --teal:#13B8A6; --teal-2:#0F9488; --cyan:#38BDF8;
  --violet:#8B7BE8; --bg:#F4F7FA; --surface:#FFFFFF; --surface-2:#EAF0F5; --border:#DCE5EC;
  --text:#172B3A; --muted:#64748B;
  --mono:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
}
*{box-sizing:border-box}
html{scroll-behavior:smooth;scroll-padding-top:4.5rem}
body{margin:0;background:var(--bg);color:var(--text);
  font:15px/1.55 Inter,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
a{color:var(--navy-2)} a:hover{color:var(--teal-2)}
code{font-family:var(--mono);font-size:.88em;background:var(--surface-2);padding:.05em .3em;border-radius:4px}
.wrap{max-width:1120px;margin:0 auto;padding:0 1.25rem}
header.top{background:var(--navy);color:#fff;position:relative;overflow:hidden}
header.top .strip{position:absolute;inset:auto 0 0 0;height:44px;opacity:.5}
header.top .inner{position:relative;padding:2.6rem 0 4rem}
.brand{display:flex;align-items:center;gap:.6rem;font-weight:600;letter-spacing:.02em;color:#BFD3E0}
.brand b{color:var(--teal)}
h1{font-size:clamp(1.6rem,3.4vw,2.35rem);line-height:1.15;margin:.9rem 0 .5rem;font-weight:650;max-width:22em}
.lede{color:#BFD3E0;max-width:46em;margin:0}
nav.bar{position:sticky;top:0;z-index:5;background:rgba(16,29,45,.96);backdrop-filter:blur(6px);
  border-bottom:1px solid #1d3550}
nav.bar .wrap{display:flex;gap:.25rem;overflow-x:auto;white-space:nowrap}
nav.bar a{color:#BFD3E0;text-decoration:none;padding:.75rem .85rem;font-size:.9rem;border-bottom:2px solid transparent}
nav.bar a:hover,nav.bar a:focus-visible{color:#fff;border-bottom-color:var(--teal);outline:none}
.notice{background:var(--surface);border:1px solid var(--border);border-left:4px solid var(--violet);
  border-radius:8px;padding:.7rem 1rem;margin:1.4rem 0 0}
section{margin:2.4rem 0}
section>h2{font-size:1.3rem;margin:0 0 .25rem;font-weight:650}
section>p.sub{color:var(--muted);margin:.1rem 0 1rem;max-width:52em}
.card{background:var(--surface);border:1px solid var(--border);border-radius:10px;padding:1.1rem 1.2rem}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:.8rem;margin-top:1.4rem}
.tile{background:var(--surface);border:1px solid var(--border);border-radius:10px;padding:.9rem 1rem;
  box-shadow:0 1px 0 rgba(16,29,45,.04)}
.tile .k{font-size:.78rem;color:var(--muted);text-transform:uppercase;letter-spacing:.06em}
.tile .v{font:600 1.6rem/1.2 var(--mono);color:var(--navy);margin-top:.2rem}
.tile .d{font-size:.82rem;color:var(--muted);margin-top:.15rem}
.controls{display:flex;flex-wrap:wrap;gap:.9rem 1.4rem;align-items:flex-end;margin-bottom:1rem}
.ctl label,.ctl .lab{display:block;font-size:.76rem;color:var(--muted);text-transform:uppercase;
  letter-spacing:.06em;margin-bottom:.25rem}
.seg{display:inline-flex;border:1px solid var(--border);border-radius:8px;overflow:hidden;background:var(--surface)}
.seg button,.pill{font:inherit;border:0;background:transparent;color:var(--text);padding:.45rem .8rem;cursor:pointer}
.seg button+button{border-left:1px solid var(--border)}
.seg button[aria-pressed=true]{background:var(--navy);color:#fff}
.seg button:hover:not([aria-pressed=true]){background:var(--surface-2)}
.pills{display:flex;gap:.4rem;flex-wrap:wrap}
.pill{border:1px solid var(--border);border-radius:999px;background:var(--surface);display:inline-flex;align-items:center;gap:.4rem}
.pill i{width:.7rem;height:.7rem;border-radius:50%;display:inline-block}
.pill[aria-pressed=true]{background:var(--navy);color:#fff;border-color:var(--navy)}
.pill:hover:not([aria-pressed=true]){background:var(--surface-2)}
button:focus-visible,select:focus-visible,a:focus-visible{outline:2px solid var(--cyan);outline-offset:2px}
select{font:inherit;padding:.42rem .6rem;border:1px solid var(--border);border-radius:8px;background:var(--surface);color:var(--text);max-width:100%}
.btn{font:inherit;padding:.45rem .85rem;border-radius:8px;border:1px solid var(--navy);background:var(--navy);color:#fff;cursor:pointer}
.btn:hover{background:var(--navy-2)}
.btn.ghost{background:var(--surface);color:var(--navy)} .btn.ghost:hover{background:var(--surface-2)}
.note{font-size:.88rem;color:var(--muted);margin:.5rem 0}
.callout{background:var(--surface-2);border-radius:8px;padding:.6rem .85rem;font-size:.9rem;margin:.6rem 0}
svg text{font-family:inherit}
.chart{width:100%;height:auto;display:block}
table{border-collapse:collapse;width:100%;font-size:.88rem}
th,td{padding:.4rem .6rem;border-bottom:1px solid var(--border);text-align:left;vertical-align:top}
th{font-size:.74rem;text-transform:uppercase;letter-spacing:.05em;color:var(--muted);background:var(--surface-2);position:sticky;top:0}
td.n{font-family:var(--mono);white-space:nowrap}
.tablewrap{overflow:auto;max-height:420px;border:1px solid var(--border);border-radius:8px;margin-top:.8rem}
.muted{color:var(--muted)}
.cols{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:1rem}
.or{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:.8rem;margin-top:.8rem}
.or .card .big{font:600 1.5rem var(--mono);color:var(--navy)}
.tag{display:inline-block;font-size:.72rem;padding:.1rem .45rem;border-radius:999px;background:var(--surface-2);color:var(--navy-2);margin-left:.4rem;vertical-align:middle}
.tag.warn{background:#EFEDFF;color:#4b3fb0}
ol.lim li{margin:.4rem 0}
.legend{display:flex;gap:1rem;flex-wrap:wrap;font-size:.85rem;margin:.5rem 0}
.legend i{display:inline-block;width:.75rem;height:.75rem;border-radius:2px;margin-right:.35rem;vertical-align:-.05rem}
footer{border-top:1px solid var(--border);padding:1.6rem 0 3rem;color:var(--muted);font-size:.85rem}
@media print{nav.bar,.controls,.btn{display:none}}
</style></head>
<body>
<header class="top">
  <svg class="strip" id="strip" aria-hidden="true" preserveAspectRatio="none"></svg>
  <div class="wrap inner">
    <div class="brand"><svg width="26" height="26" viewBox="0 0 26 26" aria-hidden="true"><path d="M5 3c0 7 16 7 16 10S5 19 5 23M21 3c0 7-16 7-16 10s16 6 16 10" fill="none" stroke="#13B8A6" stroke-width="1.8" stroke-linecap="round"/><path d="M7 8h12M6 13h14M7 18h12" stroke="#38BDF8" stroke-width="1.2" opacity=".8"/></svg><span>Viro<b>trace</b></span></div>
    <h1 id="title"></h1>
    <p class="lede">HIV-1 <em>pol</em> drug-resistance mutations and lineages in deposited sequences, one per individual, with every denominator shown. Built from public GenBank records.</p>
  </div>
</header>
<nav class="bar" aria-label="Sections"><div class="wrap">
  <a href="#resistance">Resistance</a><a href="#clusters">Clusters</a><a href="#data">Data and QC</a>
  <a href="#trees">Trees</a><a href="#limitations">Limitations</a><a href="#refs">References</a>
</div></nav>
<main class="wrap">
  <div class="tiles" id="tiles"></div>
  <p class="notice"><strong>Not a clinical tool.</strong> These are descriptive statistics from publicly deposited sequences. They describe published sequences, not the epidemic, and no result should inform the care of an individual patient.</p>
  <noscript><p class="notice">This page needs JavaScript to draw the charts. The same numbers are in <code>data.json</code>.</p></noscript>

  <section id="resistance">
    <h2>Drug-resistance prevalence</h2>
    <p class="sub">Proportion of sequences carrying at least one surveillance drug-resistance mutation (SDRM), by drug class. Bars show 95% Wilson intervals. Groups below the minimum size are suppressed.</p>
    <div class="card">
      <div class="controls">
        <div class="ctl"><span class="lab">Mutation list</span><div class="seg" id="roleSeg" role="group" aria-label="Mutation list"></div></div>
        <div class="ctl"><label for="kind">Group by</label><select id="kind"></select></div>
        <div class="ctl"><span class="lab">Drug class</span><div class="pills" id="clsPills" role="group" aria-label="Drug class"></div></div>
        <div class="ctl"><label for="sort">Order</label><select id="sort"><option value="natural">Natural</option><option value="n">Sample size</option><option value="p">Proportion</option></select></div>
        <div class="ctl"><button class="btn ghost" id="csv" type="button">Download CSV</button></div>
      </div>
      <div id="roleNote"></div>
      <div id="forest"></div>
      <div class="tablewrap"><table id="matrix"></table></div>
    </div>
    <div class="card" id="trendCard" style="margin-top:1rem"></div>
  </section>

  <section id="clusters">
    <h2>Clusters</h2>
    <p class="sub">A descriptive statistic from the tree set at a fixed TN93 distance. Clusters are not transmission chains.</p>
    <div id="clusterBody"></div>
  </section>

  <section id="data">
    <h2>Data and quality control</h2>
    <p class="sub">What went in, what was removed and why, and what the tree set looks like.</p>
    <div class="cols" id="dataBody"></div>
  </section>

  <section id="trees">
    <h2>Trees and lineages</h2>
    <p class="sub">Temporal signal decides whether a tree is shown dated or as a divergence tree.</p>
    <div id="treeBody"></div>
  </section>

  <section id="limitations">
    <h2>Limitations</h2>
    <p class="sub">Read these before quoting a number.</p>
    <div class="card"><ol class="lim" id="limList"></ol></div>
  </section>

  <section id="refs">
    <h2>References</h2>
    <div class="card"><ol id="refList"></ol></div>
  </section>
</main>
<footer><div class="wrap" id="foot"></div></footer>

<script id="virotrace-data" type="application/json">__DATA__</script>
<script>
(function(){
"use strict";
var D=JSON.parse(document.getElementById('virotrace-data').textContent);
var $=function(s,r){return (r||document).querySelector(s)};
var esc=function(s){return String(s==null?'':s).replace(/[&<>"']/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]})};
var pct=function(v,d){return v==null?'':(v*100).toFixed(d==null?1:d)+'%'};
var fmt=function(n){return n==null?'':Number(n).toLocaleString('en-GB')};
var COLOR={NRTI:'#13B8A6',NNRTI:'#38BDF8',PI:'#8B7BE8',INSTI:'#173B4D'};
var ROLE_LABEL={headline:'WHO 2009 (headline)',sensitivity:'Stanford (all positions)',sensitivity_headline_denominator:'Stanford (same people)'};
var ROLE_ORDER=['headline','sensitivity','sensitivity_headline_denominator'];
var scope=D.scope.charAt(0).toUpperCase()+D.scope.slice(1);
document.title=document.title||'Virotrace';
$('#title').textContent=scope+' HIV-1 drug resistance and lineage explorer';

/* decorative alignment strip: a deterministic pseudo-random run of base colours */
(function(){var s=$('#strip'),w=1200,h=44,x=0,seed=7,out='';var cols=['#13B8A6','#38BDF8','#8B7BE8','#173B4D'];
  s.setAttribute('viewBox','0 0 '+w+' '+h);
  while(x<w){seed=(seed*1103515245+12345)&0x7fffffff;var c=cols[seed%4],bh=6+(seed>>4)%30;
    out+='<rect x="'+x+'" y="'+(h-bh)+'" width="6" height="'+bh+'" rx="1" fill="'+c+'" opacity="'+(0.25+((seed>>9)%5)/10)+'"/>';x+=8}
  s.innerHTML=out})();

/* ---------- summary tiles ---------- */
(function(){
  var t=[],qc=D.qc,ts=D.tree_set,cl=D.clusters;
  if(qc&&qc.n_passed!=null)t.push(['Sequences passing QC',fmt(qc.n_passed),'before population policies']);
  var after=qc&&qc.corpus_after_qc_and_clone_policy;
  if(after&&after.n_total!=null)t.push(['Individuals analysed',fmt(after.n_total),'one sequence per person']);
  if(ts){var wf=ts.waterfall||{};var n=wf.after_neighbour_screen!=null?wf.after_neighbour_screen:wf.retained;
    if(n!=null)t.push(['Tree set',fmt(n),'strictly dated individuals'])}
  if(cl&&cl.proportion_clustered!=null)t.push(['In a cluster',pct(cl.proportion_clustered),fmt(cl.n_clusters)+' clusters']);
  var head=D.prevalence.filter(function(r){return r.list_role==='headline'&&r.stratum_kind==='overall'&&r.status==='estimated'});
  if(head.length){var m=head.reduce(function(a,b){return b.denominator>a.denominator?b:a});
    t.push(['Largest overall denominator',fmt(m.denominator),m.drug_class+' on the WHO 2009 list'])}
  $('#tiles').innerHTML=t.map(function(x){return '<div class="tile"><div class="k">'+esc(x[0])+'</div><div class="v">'+esc(x[1])+'</div><div class="d">'+esc(x[2])+'</div></div>'}).join('');
})();

/* ---------- resistance explorer ---------- */
var roles=ROLE_ORDER.filter(function(r){return D.prevalence.some(function(x){return x.list_role===r})});
var classes=['NRTI','NNRTI','PI','INSTI'].filter(function(c){return D.prevalence.some(function(x){return x.drug_class===c})});
var st={role:roles[0],kind:null,cls:classes[0],sort:null};
function kindsFor(role){var k=[];D.prevalence.forEach(function(r){if(r.list_role===role&&k.indexOf(r.stratum_kind)<0)k.push(r.stratum_kind)});
  var first=['overall'];return k.sort(function(a,b){return (first.indexOf(b)>=0)-(first.indexOf(a)>=0)})}
function kindLabel(k){return k==='overall'?'Overall':k==='naive_year_bin'?'Treatment-naive studies, by year bin':'By '+k.replace(/_/g,' ')}
var YEARISH=function(k){return /year/.test(k)};

function setRole(role){st.role=role;var ks=kindsFor(role);if(ks.indexOf(st.kind)<0)st.kind=ks[0];st.sort=null;drawControls();drawExplorer()}
function drawControls(){
  $('#roleSeg').innerHTML=roles.map(function(r){return '<button type="button" data-role="'+r+'" aria-pressed="'+(r===st.role)+'">'+esc(ROLE_LABEL[r]||r)+'</button>'}).join('');
  var ks=kindsFor(st.role);
  $('#kind').innerHTML=ks.map(function(k){return '<option value="'+esc(k)+'"'+(k===st.kind?' selected':'')+'>'+esc(kindLabel(k))+'</option>'}).join('');
  $('#clsPills').innerHTML=classes.map(function(c){return '<button type="button" class="pill" data-cls="'+c+'" aria-pressed="'+(c===st.cls)+'"><i style="background:'+COLOR[c]+'"></i>'+c+'</button>'}).join('');
  $('#sort').value=st.sort||(YEARISH(st.kind)||st.kind==='overall'?'natural':'n');
}
$('#roleSeg').addEventListener('click',function(e){var b=e.target.closest('button');if(b)setRole(b.dataset.role)});
$('#clsPills').addEventListener('click',function(e){var b=e.target.closest('button');if(b){st.cls=b.dataset.cls;drawControls();drawExplorer()}});
$('#kind').addEventListener('change',function(e){st.kind=e.target.value;st.sort=null;drawControls();drawExplorer()});
$('#sort').addEventListener('change',function(e){st.sort=e.target.value;drawExplorer()});
$('#csv').addEventListener('click',function(){
  var rows=D.prevalence.filter(function(r){return r.list_role===st.role&&r.stratum_kind===st.kind});
  var cols=['list_role','list_version','stratum_kind','stratum','drug_class','numerator','denominator','proportion','ci_low','ci_high','status'];
  var q=function(v){v=v==null?'':String(v);return /[",\n]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v};
  var text=[cols.join(',')].concat(rows.map(function(r){return cols.map(function(c){return q(r[c])}).join(',')})).join('\n');
  var a=document.createElement('a');a.href=URL.createObjectURL(new Blob([text],{type:'text/csv'}));
  a.download='virotrace_'+D.scope+'_'+st.role+'_'+st.kind+'.csv';document.body.appendChild(a);a.click();a.remove()});

function roleNote(){
  var v=D.list_versions[st.role]||'';var n=[];
  n.push('<div class="callout"><strong>'+esc(ROLE_LABEL[st.role]||st.role)+'</strong>'+(v?' <span class="muted">version '+esc(v)+'</span>':'')+'. ');
  if(st.role==='headline')n.push('The WHO 2009 surveillance list is the headline. Every sequence is counted at the positions it covers, so figures are lower bounds.');
  if(st.role==='sensitivity')n.push('Only sequences covering every position of this list are counted, so denominators are smaller than the headline and the percentages are not comparable to it.');
  if(st.role==='sensitivity_headline_denominator')n.push('The same people as the headline with the Stanford mutations called at the positions each sequence covers. This is the like-for-like comparison with the headline.');
  n.push('</div>');
  if(st.kind==='naive_year_bin')n.push('<div class="callout">These are pooled treatment-naive <em>studies</em>, not a transmitted-resistance trend: study composition changes between bins, and the earliest bins include natural polymorphisms and small denominators.</div>');
  if(st.kind==='overall'&&D.scope==='africa')n.push('<div class="callout">The overall row pools every country by sequence count. It is dominated by the best-sampled countries and is not an estimate for the region.</div>');
  $('#roleNote').innerHTML=n.join('');
}

function drawExplorer(){
  roleNote();
  var rows=D.prevalence.filter(function(r){return r.list_role===st.role&&r.stratum_kind===st.kind});
  var strata=[];rows.forEach(function(r){if(strata.indexOf(r.stratum)<0)strata.push(r.stratum)});
  var items=strata.map(function(s){return {s:s,r:rows.find(function(x){return x.stratum===s&&x.drug_class===st.cls})}}).filter(function(x){return x.r});
  var sort=$('#sort').value;
  if(sort==='n')items.sort(function(a,b){return (b.r.denominator||0)-(a.r.denominator||0)});
  else if(sort==='p')items.sort(function(a,b){return (b.r.proportion==null?-1:b.r.proportion)-(a.r.proportion==null?-1:a.r.proportion)});
  else items.sort(function(a,b){return a.s<b.s?-1:a.s>b.s?1:0});
  var cut=items.length>40;if(cut)items=items.slice().sort(function(a,b){return (b.r.denominator||0)-(a.r.denominator||0)}).slice(0,40);
  var est=items.filter(function(x){return x.r.status==='estimated'});
  var max=Math.max.apply(null,[0.1].concat(est.map(function(x){return x.r.ci_high||0})));
  max=Math.min(1,Math.ceil(max*20)/20);
  var W=900,L=230,R=290,rowH=28,top=26,H=top+items.length*rowH+34,pw=W-L-R;
  var X=function(v){return L+pw*v/max};
  var svg='<svg class="chart" viewBox="0 0 '+W+' '+H+'" role="img" aria-label="Prevalence of '+esc(st.cls)+' surveillance mutations by group, with 95% intervals">';
  var step=max>0.5?0.2:max>0.25?0.1:0.05;
  for(var t=0;t<=max+1e-9;t+=step){svg+='<line x1="'+X(t)+'" x2="'+X(t)+'" y1="'+(top-6)+'" y2="'+(H-28)+'" stroke="#DCE5EC"/><text x="'+X(t)+'" y="'+(H-12)+'" text-anchor="middle" font-size="11" fill="#64748B">'+Math.round(t*100)+'%</text>'}
  items.forEach(function(it,i){var y=top+i*rowH+rowH/2,r=it.r,lab=it.s.length>30?it.s.slice(0,29)+'\u2026':it.s;
    svg+='<text x="'+(L-10)+'" y="'+(y+4)+'" text-anchor="end" font-size="12.5" fill="#172B3A"><title>'+esc(it.s)+'</title>'+esc(lab)+'</text>';
    if(r.status==='estimated'){
      svg+='<line x1="'+X(r.ci_low)+'" x2="'+X(r.ci_high)+'" y1="'+y+'" y2="'+y+'" stroke="'+COLOR[st.cls]+'" stroke-width="3" stroke-linecap="round" opacity=".55"/>'+
        '<circle cx="'+X(r.proportion)+'" cy="'+y+'" r="5" fill="'+COLOR[st.cls]+'" stroke="#101D2D" stroke-width="1"><title>'+esc(it.s)+': '+pct(r.proportion)+' ('+pct(r.ci_low)+' to '+pct(r.ci_high)+'), '+r.numerator+'/'+r.denominator+'</title></circle>'+
        '<text x="'+(W-R+14)+'" y="'+(y+4)+'" font-size="12" fill="#172B3A" font-family="ui-monospace,Menlo,monospace">'+pct(r.proportion)+' <tspan fill="#64748B">('+pct(r.ci_low)+'\u2013'+pct(r.ci_high)+')  '+fmt(r.numerator)+'/'+fmt(r.denominator)+'</tspan></text>';
    }else{
      svg+='<text x="'+L+'" y="'+(y+4)+'" font-size="12" fill="#64748B" font-style="italic">'+(r.status==='n_too_small'?'suppressed, n='+fmt(r.denominator):'no list positions')+'</text>'}
  });
  svg+='</svg>';
  $('#forest').innerHTML=(items.length?svg:'<p class="note">No estimates for this class in this grouping.</p>')+(cut?'<p class="note">Showing the 40 largest groups; the table below and the CSV have all of them.</p>':'');
  /* matrix: every class, shaded by proportion */
  var head='<thead><tr><th>Group</th>'+classes.map(function(c){return '<th><span style="color:'+COLOR[c]+'">\u25A0</span> '+c+'</th>'}).join('')+'</tr></thead>';
  var all=strata.slice();
  if(sort==='n')all.sort(function(a,b){return den(b)-den(a)});else if(sort==='p')all.sort(function(a,b){return pr(b)-pr(a)});else all.sort();
  function cell(s,c){return rows.find(function(x){return x.stratum===s&&x.drug_class===c})}
  function den(s){var r=cell(s,st.cls);return r&&r.denominator||0}
  function pr(s){var r=cell(s,st.cls);return r&&r.proportion!=null?r.proportion:-1}
  var body='<tbody>'+all.map(function(s){return '<tr><td>'+esc(s)+'</td>'+classes.map(function(c){var r=cell(s,c);
    if(!r)return '<td class="muted">not assessed</td>';
    if(r.status!=='estimated')return '<td class="muted">'+(r.status==='n_too_small'?'n='+fmt(r.denominator)+' (suppressed)':'no list positions')+'</td>';
    var a=Math.min(0.55,r.proportion*1.6);
    return '<td class="n" style="background:rgba(19,184,166,'+a.toFixed(2)+')">'+pct(r.proportion)+' <span class="muted">'+fmt(r.numerator)+'/'+fmt(r.denominator)+'</span></td>'}).join('')+'</tr>'}).join('')+'</tbody>';
  $('#matrix').innerHTML=head+body;
  drawTrends();
}
function pText(p){return p==null?'':p===0?'< 1e-300':String(p)}
function drawTrends(){
  var t=D.trends||[];var c=$('#trendCard');
  if(!t.length||st.role!=='headline'){c.style.display='none';return}
  c.style.display='';
  c.innerHTML='<h3 style="margin:.1rem 0 .3rem;font-size:1rem">Trend tests (headline list)</h3><p class="note">A trend can reflect changing study composition rather than the epidemic. Very small p-values come from very large pooled samples and are not a measure of importance.</p><div class="tablewrap" style="max-height:none"><table><thead><tr><th>Class</th><th>Bins</th><th>First</th><th>Last</th><th>z</th><th>p (two-sided)</th><th>Note</th></tr></thead><tbody>'+
    t.map(function(r){return '<tr><td>'+esc(r.drug_class)+'</td><td class="n">'+esc(r.n_bins)+'</td><td>'+esc(r.first_bin)+'</td><td>'+esc(r.last_bin)+'</td><td class="n">'+esc(r.z)+'</td><td class="n">'+esc(pText(r.p_two_sided))+'</td><td>'+esc(r.note)+'</td></tr>'}).join('')+'</tbody></table></div>';
}

/* ---------- clusters ---------- */
function orText(b){return b&&b.odds_ratio!=null?{v:b.odds_ratio,ci:b.ci_low+' to '+b.ci_high}:null}
(function(){
  var cl=D.clusters,el=$('#clusterBody');
  if(!cl){el.innerHTML='<p class="note">Not generated for this build.</p>';return}
  var s=cl.sdrm||{},c=s.clustered||{},u=s.unclustered||{};
  var out='<div class="card"><p style="margin-top:0">At a TN93 distance of <strong>'+esc(cl.threshold)+'</strong>, <strong>'+fmt(cl.n_clustered)+'</strong> of '+fmt(cl.n_sequences)+' tree-set sequences ('+pct(cl.proportion_clustered)+') fall in <strong>'+fmt(cl.n_clusters)+'</strong> clusters; the largest has '+fmt(cl.largest_cluster)+'.</p>';
  var bars=[['In a cluster',c],['Not in a cluster',u]];
  out+='<svg class="chart" viewBox="0 0 900 96" role="img" aria-label="Proportion of sequences with an SDRM, in clusters and not">';
  bars.forEach(function(b,i){var y=14+i*40,p=b[1].proportion||0;
    out+='<text x="150" y="'+(y+15)+'" text-anchor="end" font-size="13" fill="#172B3A">'+b[0]+'</text><rect x="162" y="'+y+'" width="560" height="22" rx="4" fill="#EAF0F5"/><rect x="162" y="'+y+'" width="'+(560*p)+'" height="22" rx="4" fill="'+(i?'#8B7BE8':'#13B8A6')+'"/><text x="732" y="'+(y+16)+'" font-size="12.5" fill="#172B3A" font-family="ui-monospace,Menlo,monospace">'+pct(b[1].proportion)+'  '+fmt(b[1].with_sdrm)+'/'+fmt(b[1].n)+'</text>'});
  out+='</svg><p class="note" style="margin-bottom:0">Share of sequences carrying an SDRM. Counts come from the tree set only.</p></div><div class="or">';
  var a=cl.sdrm_cluster_odds_ratio_adjusted_for_study,ao=orText(a);
  if(ao)out+='<div class="card"><div class="k muted">Study-adjusted odds ratio <span class="tag">headline</span></div><div class="big">'+esc(ao.v)+'</div><div class="muted">95% interval '+esc(ao.ci)+'</div><p class="note">Mantel-Haenszel with each study a stratum; '+esc(a.n_studies_informative)+' of '+esc(a.n_studies)+' studies contribute. An interval that includes 1 means no detectable difference, not evidence of equal transmission.</p></div>';
  var tc=orText(cl.sdrm_cluster_odds_ratio_adjusted_for_treatment_context);
  if(tc)out+='<div class="card"><div class="k muted">Treatment-context-adjusted <span class="tag warn">confounded</span></div><div class="big">'+esc(tc.v)+'</div><div class="muted">95% interval '+esc(tc.ci)+'</div><p class="note">Confounded by study structure: a few densely sampled cohorts with little resistance dominate the clusters, and many sequences have no treatment annotation. Compare it with the study-adjusted figure; it is not a finding.</p></div>';
  var alt=cl.sensitivity_excluding_studies,ar=alt&&orText(alt.sdrm_cluster_odds_ratio_adjusted_for_study);
  if(ar)out+='<div class="card"><div class="k muted">Study-adjusted, studies removed <span class="tag">sensitivity</span></div><div class="big">'+esc(ar.v)+'</div><div class="muted">95% interval '+esc(ar.ci)+'</div><p class="note">With '+esc((alt.excluded_studies||[]).join(', '))+' removed and the rest re-clustered ('+fmt(alt.n_sequences_removed)+' sequences).</p></div>';
  el.innerHTML=out+'</div>';
})();

/* ---------- data and QC ---------- */
function hbars(title,pairs,color,note){
  if(!pairs.length)return '';
  var max=Math.max.apply(null,pairs.map(function(p){return p[1]}))||1,rowH=24,H=pairs.length*rowH+8;
  var s='<div class="card"><h3 style="margin:.1rem 0 .5rem;font-size:1rem">'+esc(title)+'</h3><svg class="chart" viewBox="0 0 520 '+H+'" role="img" aria-label="'+esc(title)+'">';
  pairs.forEach(function(p,i){var y=4+i*rowH,w=Math.max(1,300*p[1]/max);
    s+='<text x="150" y="'+(y+15)+'" text-anchor="end" font-size="12" fill="#172B3A">'+esc(String(p[0]).slice(0,26))+'</text><rect x="158" y="'+(y+3)+'" width="'+w+'" height="16" rx="3" fill="'+color+'"/><text x="'+(166+w)+'" y="'+(y+15)+'" font-size="11.5" fill="#64748B" font-family="ui-monospace,Menlo,monospace">'+fmt(p[1])+'</text>'});
  return s+'</svg>'+(note?'<p class="note" style="margin-bottom:0">'+note+'</p>':'')+'</div>';
}
function columns(title,pairs,color){
  if(!pairs.length)return '';
  var max=Math.max.apply(null,pairs.map(function(p){return p[1]}))||1,W=520,H=170,bw=Math.min(28,(W-40)/pairs.length-3);
  var s='<div class="card"><h3 style="margin:.1rem 0 .5rem;font-size:1rem">'+esc(title)+'</h3><svg class="chart" viewBox="0 0 '+W+' '+H+'" role="img" aria-label="'+esc(title)+'">';
  pairs.forEach(function(p,i){var x=30+i*(bw+3),h=(H-40)*p[1]/max;
    s+='<rect x="'+x+'" y="'+(H-24-h)+'" width="'+bw+'" height="'+h+'" rx="2" fill="'+color+'"><title>'+esc(p[0])+': '+fmt(p[1])+'</title></rect>';
    if(pairs.length<=16||i%2===0)s+='<text x="'+(x+bw/2)+'" y="'+(H-8)+'" text-anchor="middle" font-size="10" fill="#64748B">'+esc(String(p[0]).slice(-2)==='00'?p[0]:String(p[0]).slice(2))+'</text>'});
  return s+'</svg></div>';
}
(function(){
  var qc=D.qc,ts=D.tree_set,out='';
  if(qc){var w=qc.waterfall||{};
    var pairs=Object.keys(w).sort().map(function(k){return [k,w[k]]});
    var pol=qc.n_removed_by_population_policy||{};Object.keys(pol).forEach(function(k){pairs.push([k.replace(/_/g,' '),pol[k]])});
    out+=hbars('Sequences removed at each QC step',pairs,'#173B4D','Clone policy <code>'+esc(qc.clone_policy)+'</code>, individual policy <code>'+esc(qc.individual_policy)+'</code>.')}
  if(ts){var sub=Object.keys(ts.by_subtype||{}).map(function(k){return [k,ts.by_subtype[k]]}).sort(function(a,b){return b[1]-a[1]});
    out+=hbars('Tree set by subtype',sub,'#13B8A6');
    var yrs=Object.keys(ts.by_year||{}).sort().map(function(k){return [k,ts.by_year[k]]});
    out+=columns('Tree set by collection year',yrs,'#38BDF8');
    var ns=ts.neighbour_screen&&ts.neighbour_screen.removed||[];
    if(ns.length)out+='<div class="card"><h3 style="margin:.1rem 0 .4rem;font-size:1rem">Neighbour screen</h3><p class="note" style="margin:0">'+ns.length+' sequences were left out of the tree because no other sequence is within the identity floor of '+pct(ts.neighbour_screen.threshold,0)+' over shared sites. They remain in the resistance prevalence.</p></div>'}
  $('#dataBody').innerHTML=out||'<p class="note">Not generated for this build.</p>';
})();

/* ---------- trees ---------- */
(function(){
  var el=$('#treeBody'),out='';
  if(D.gates&&D.gates.length){
    out+='<div class="card"><h3 style="margin:.1rem 0 .3rem;font-size:1rem">Temporal signal</h3><div class="tablewrap" style="max-height:none"><table><thead><tr><th>Tree</th><th>Tips</th><th>Slope</th><th>R\u00B2</th><th>Permutation p</th><th>Gate</th><th>Shown as</th></tr></thead><tbody>'+
      D.gates.map(function(g){return '<tr><td>'+esc(g.label)+'</td><td class="n">'+fmt(g.n)+'</td><td class="n">'+esc(g.slope==null?'':Number(g.slope).toPrecision(3))+'</td><td class="n">'+esc(g.r_squared==null?'':Number(g.r_squared).toFixed(3))+'</td><td class="n">'+esc(g.permutation_p==null?'':g.permutation_p)+'</td><td>'+(g.passed?'pass':'fail')+'</td><td>'+esc(g.result)+'</td></tr>'}).join('')+'</tbody></table></div><p class="note" style="margin-bottom:0">Root-to-tip regression at the best root with a date-randomisation test. A tree that fails is shown as a midpoint-rooted divergence tree.</p></div>'}
  var st2=D.subtype_tree;
  if(st2){var c=st2.concordance||{},keys=[['agree','#13B8A6'],['disagree','#8B7BE8'],['not_comparable','#DCE5EC'],['unresolved','#38BDF8']];
    var tot=keys.reduce(function(a,k){return a+(c[k[0]]||0)},0)||1,x=0,bar='';
    keys.forEach(function(k){var w=900*(c[k[0]]||0)/tot;bar+='<rect x="'+x+'" y="0" width="'+w+'" height="26" fill="'+k[1]+'"><title>'+k[0].replace('_',' ')+': '+fmt(c[k[0]]||0)+'</title></rect>';x+=w});
    out+='<div class="card" style="margin-top:1rem"><h3 style="margin:.1rem 0 .3rem;font-size:1rem">Second subtype call (tree-based)</h3><p class="note" style="margin-top:0">'+fmt(st2.n_queries)+' tree-set sequences placed in a tree of '+fmt(st2.n_references)+' references. '+(st2.support_used?'Calls needed bootstrap support on the joining branch.':'This tree has no bootstrap support values, so no support threshold was applied and calls are weaker than a support-filtered check.')+'</p><svg class="chart" viewBox="0 0 900 26" role="img" aria-label="Agreement of the tree-based subtype call with the distance-based call">'+bar+'</svg><div class="legend">'+keys.map(function(k){return '<span><i style="background:'+k[1]+'"></i>'+k[0].replace('_',' ')+' '+fmt(c[k[0]]||0)+'</span>'}).join('')+'</div><p class="note" style="margin-bottom:0">Unresolved is common because few references represent each subtype, so this check can confirm a call but not rule one out. Disagreements are flags, not corrections.</p></div>'}
  if(D.trees&&D.trees.length){
    out+='<div class="card" style="margin-top:1rem"><h3 style="margin:.1rem 0 .3rem;font-size:1rem">Interactive trees (Auspice)</h3><p class="note" style="margin-top:0">Download a file and drop it on <a href="https://auspice.us" rel="noopener">auspice.us</a>; the file is read in your browser and not uploaded. Colourings include resistance by class, subtype, country and treatment context.</p><div class="tablewrap" style="max-height:none"><table><thead><tr><th>Tree</th><th>File</th><th>Size</th><th></th></tr></thead><tbody>'+
      D.trees.map(function(t){return '<tr><td>'+esc(t.label)+'</td><td><a href="'+esc(t.file)+'" download>'+esc(t.name)+'</a></td><td class="n">'+(t.bytes/1048576).toFixed(1)+' MB</td><td>'+(t.nextstrain?'<a href="'+esc(t.nextstrain)+'" rel="noopener">Open on nextstrain.org</a>':'')+'</td></tr>'}).join('')+'</tbody></table></div></div>'}
  el.innerHTML=out||'<p class="note">Not generated for this build.</p>';
})();

/* ---------- limitations, references, footer ---------- */
$('#limList').innerHTML=D.limitations.map(function(t){return '<li>'+esc(t)+'</li>'}).join('');
$('#refList').innerHTML=D.references.map(function(t){return '<li>'+esc(t)+'</li>'}).join('');
var vs=Object.keys(D.list_versions).map(function(k){return (ROLE_LABEL[k]||k)+' '+D.list_versions[k]}).join('; ');
$('#foot').innerHTML='Virotrace explorer, scope <code>'+esc(D.scope)+'</code>, generated '+esc(D.generated)+'. Mutation lists: '+esc(vs)+'. Minimum group size '+esc(D.min_stratum_n)+'. Aggregate tables only: no accessions or per-sequence rows are in this page. The same data is in <code>data.json</code>.';

setRole(st.role);
})();
</script>
</body></html>
"""
