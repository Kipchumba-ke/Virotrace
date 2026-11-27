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
#mapCard .maprow{display:grid;grid-template-columns:minmax(260px,560px) 1fr;gap:1.2rem;align-items:start}
@media (max-width:760px){#mapCard .maprow{grid-template-columns:1fr}}
.mapsvg{width:100%;height:auto;max-height:640px;display:block}
.mapsvg path{stroke:#fff;stroke-width:.8;cursor:pointer;transition:opacity .1s}
.mapsvg path:hover{opacity:.82} .mapsvg path:focus-visible{outline:none;stroke:#101D2D;stroke-width:2.2}
.mapsvg path.sel{stroke:#101D2D;stroke-width:2.4}
.ramp{height:10px;border-radius:5px;border:1px solid var(--border)}
.swatch{display:inline-block;width:.8rem;height:.8rem;border-radius:2px;vertical-align:-.1rem;margin-right:.35rem;border:1px solid var(--border)}
.tvbar{display:flex;flex-wrap:wrap;gap:.6rem .9rem;align-items:flex-end;margin-bottom:.7rem}
.tvbar input[type=search]{font:inherit;padding:.42rem .6rem;border:1px solid var(--border);border-radius:8px;width:11rem}
.tvbar .btn{padding:.4rem .7rem}
.tvwrap{position:relative;border:1px solid var(--border);border-radius:8px;background:#fff;overflow:hidden}
#tvCanvas{display:block;width:100%;height:640px;cursor:grab;touch-action:none}
#tvCanvas.drag{cursor:grabbing}
#tvTip{position:absolute;pointer-events:none;background:var(--navy);color:#fff;border-radius:6px;padding:.4rem .6rem;font-size:.8rem;line-height:1.35;max-width:300px;display:none;z-index:2}
#tvTip b{font-family:var(--mono);font-weight:600}
.legend button{font:inherit;border:1px solid transparent;background:transparent;color:inherit;border-radius:999px;padding:.1rem .5rem;cursor:pointer}
.legend button:hover{background:var(--surface-2)} .legend button[aria-pressed=true]{border-color:var(--navy);background:var(--surface-2)}
#tvPin table{margin-top:.5rem}
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
    <div class="card" id="mapCard" style="margin-top:1rem"></div>
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
  var an=D.annotation;
  if(an&&an.n_sequences)t.push(['Treatment context known',pct(an.n_annotated/an.n_sequences),fmt(an.n_studies_annotated)+' of '+fmt(an.n_studies)+' studies annotated']);
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
  drawTrends();drawMap();
}

/* ---------- country map ---------- */
function lerp(a,b,t){var h=function(x,i){return parseInt(x.slice(1+2*i,3+2*i),16)};
  return 'rgb('+[0,1,2].map(function(i){return Math.round(h(a,i)+(h(b,i)-h(a,i))*t)}).join(',')+')'}
var mapState={metric:'p',sel:null};
function drawMap(){
  var el=$('#mapCard'),M=D.map;
  if(!M){el.style.display='none';return}
  el.style.display='';
  var rows=D.prevalence.filter(function(r){return r.list_role===st.role&&r.stratum_kind==='country'&&r.drug_class===st.cls});
  var byId={};rows.forEach(function(r){var id=D.country_ids[r.stratum];if(id)byId[id]=r});
  var est=Object.keys(byId).map(function(k){return byId[k]}).filter(function(r){return r.status==='estimated'});
  var maxP=Math.max.apply(null,[0.05].concat(est.map(function(r){return r.proportion}))),
      maxN=Math.max.apply(null,[1].concat(Object.keys(byId).map(function(k){return byId[k].denominator||0})));
  var isP=mapState.metric==='p';
  function fill(f){var r=byId[f.id];if(!r)return '#EAF0F5';
    if(isP){if(r.status!=='estimated')return 'url(#hatch)';return lerp('#E3F6F3','#0A6A60',Math.min(1,r.proportion/maxP))}
    return lerp('#E4ECF2','#173B4D',Math.log(1+(r.denominator||0))/Math.log(1+maxN))}
  function label(f){var r=byId[f.id];if(!r)return f.name+': no sequences in this grouping';
    if(r.status!=='estimated')return f.name+': suppressed, n='+fmt(r.denominator);
    return f.name+': '+pct(r.proportion)+' ('+pct(r.ci_low)+' to '+pct(r.ci_high)+'), '+fmt(r.numerator)+'/'+fmt(r.denominator)}
  var svg='<svg class="mapsvg" viewBox="0 0 '+M.width+' '+M.height+'" role="group" aria-label="Map of Africa shaded by '+(isP?esc(st.cls)+' prevalence':'sequence count')+'"><defs><pattern id="hatch" width="7" height="7" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="7" height="7" fill="#EFEDFF"/><line x1="0" y1="0" x2="0" y2="7" stroke="#8B7BE8" stroke-width="2"/></pattern></defs>';
  M.features.forEach(function(f){
    svg+='<path d="'+f.d+'" data-id="'+esc(f.id)+'" fill="'+fill(f)+'" class="'+(mapState.sel===f.id?'sel':'')+'" tabindex="0" role="button" aria-label="'+esc(label(f))+'"><title>'+esc(label(f))+'</title></path>'});
  svg+='</svg>';
  var legend=isP
    ?'<div class="ramp" style="background:linear-gradient(90deg,#E3F6F3,#0A6A60)"></div><div style="display:flex;justify-content:space-between;font-size:.8rem" class="muted"><span>0%</span><span>'+pct(maxP,0)+'</span></div><p class="note"><span class="swatch" style="background:repeating-linear-gradient(45deg,#8B7BE8 0 2px,#EFEDFF 2px 5px)"></span>suppressed (small group) &nbsp; <span class="swatch" style="background:#EAF0F5"></span>no sequences</p>'
    :'<div class="ramp" style="background:linear-gradient(90deg,#E4ECF2,#173B4D)"></div><div style="display:flex;justify-content:space-between;font-size:.8rem" class="muted"><span>1</span><span>'+fmt(maxN)+' (log scale)</span></div><p class="note"><span class="swatch" style="background:#EAF0F5"></span>no sequences</p>';
  var sel=mapState.sel&&M.features.filter(function(f){return f.id===mapState.sel})[0],detail='<p class="note">Click a country for its estimates in every class.</p>';
  if(sel){var nm=sel.name;var allr=D.prevalence.filter(function(r){return r.list_role===st.role&&r.stratum_kind==='country'&&D.country_ids[r.stratum]===sel.id});
    detail='<h3 style="margin:.2rem 0;font-size:1rem">'+esc(nm)+'</h3>'+(allr.length?'<table><tbody>'+classes.map(function(c){var r=allr.filter(function(x){return x.drug_class===c})[0];
      return '<tr><td><span class="swatch" style="background:'+COLOR[c]+'"></span>'+c+'</td><td class="n">'+(!r?'not assessed':r.status!=='estimated'?'suppressed (n='+fmt(r.denominator)+')':pct(r.proportion)+' <span class="muted">('+pct(r.ci_low)+'\u2013'+pct(r.ci_high)+') '+fmt(r.numerator)+'/'+fmt(r.denominator)+'</span>')+'</td></tr>'}).join('')+'</tbody></table>':'<p class="note">No sequences from this country in the corpus.</p>')}
  var miss=D.unmapped_countries&&D.unmapped_countries.length?'<p class="note">Not on the map (no outline for the name): '+esc(D.unmapped_countries.join(', '))+'. They are in the table above.</p>':'';
  el.innerHTML='<div style="display:flex;flex-wrap:wrap;gap:.6rem 1.2rem;justify-content:space-between;align-items:center"><h3 style="margin:.1rem 0;font-size:1rem">Map by country</h3><div class="seg" id="mapSeg" role="group" aria-label="Map shading"><button type="button" data-m="p" aria-pressed="'+isP+'">Prevalence</button><button type="button" data-m="n" aria-pressed="'+!isP+'">Sequences</button></div></div>'+
    '<p class="note" style="margin-top:.3rem">Uses the mutation list and drug class chosen above ('+esc(ROLE_LABEL[st.role]||st.role)+', '+esc(st.cls)+'). Shading is the proportion with an SDRM; read it beside the sequence count, because a dark country with few sequences is an uncertain estimate. Prevalence describes deposited sequences, not the epidemic.</p>'+
    '<div class="maprow"><div>'+svg+'</div><div>'+legend+detail+miss+'</div></div>';
  $('#mapSeg').addEventListener('click',function(e){var b=e.target.closest('button');if(b){mapState.metric=b.dataset.m;drawMap()}});
  var pick=function(e){var p=e.target.closest('path');if(!p)return;if(e.type==='keydown'&&e.key!=='Enter'&&e.key!==' ')return;e.preventDefault();mapState.sel=p.dataset.id;drawMap();var q=$('#mapCard path[data-id="'+p.dataset.id+'"]');if(q)q.focus()};
  el.querySelector('svg').addEventListener('click',pick);el.querySelector('svg').addEventListener('keydown',pick);
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
  var an=D.annotation;
  if(an){var order=[['naive','#13B8A6'],['experienced','#173B4D'],['mixed','#8B7BE8'],['unknown','#38BDF8'],['unannotated','#DCE5EC']],bc=an.by_context||{},x=0,bar='';
    order.forEach(function(o){var w=520*(bc[o[0]]||0)/an.n_sequences;bar+='<rect x="'+x+'" y="0" width="'+w+'" height="24" fill="'+o[1]+'"><title>'+o[0]+': '+fmt(bc[o[0]]||0)+'</title></rect>';x+=w});
    out+='<div class="card"><h3 style="margin:.1rem 0 .5rem;font-size:1rem">Treatment-context annotation</h3><p class="note" style="margin-top:0"><strong>'+pct(an.n_annotated/an.n_sequences)+'</strong> of '+fmt(an.n_sequences)+' QC-passed sequences ('+fmt(an.n_studies_annotated)+' of '+fmt(an.n_studies)+' studies) have a hand-curated treatment context.</p><svg class="chart" viewBox="0 0 520 24" role="img" aria-label="Sequences by study treatment context">'+bar+'</svg><div class="legend">'+order.map(function(o){return '<span><i style="background:'+o[1]+'"></i>'+o[0]+' '+fmt(bc[o[0]]||0)+'</span>'}).join('')+'</div><p class="note" style="margin-bottom:0">Contexts come from reading each study, never from mutations. Mixed and unknown studies are annotated but not placed in the naive or experienced strata; unannotated sequences are left out of treatment strata.</p></div>'}
  $('#dataBody').innerHTML=out||'<p class="note">Not generated for this build.</p>';
})();

/* ---------- trees ---------- */
(function(){
  var el=$('#treeBody'),out='<div class="card" id="tv" style="margin-bottom:1rem"></div>';
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
  el.innerHTML=out;
})();


/* ---------- tree viewer: reads an Auspice v2 JSON, draws it on a canvas ---------- */
var TV_PALETTE=['#13B8A6','#8B7BE8','#38BDF8','#F59E0B','#E5566D','#173B4D','#7ACBBF','#B7ABF2','#8ED3F7','#F7C46A','#EE9AA7','#5C7A8C'];
var TV_FIXED={yes:'#8B7BE8',no:'#AEBBC8',unknown:'#CBD5E1',unannotated:'#CBD5E1',unclustered:'#CBD5E1',none:'#CBD5E1',unassigned:'#CBD5E1','?':'#CBD5E1'};
var TV_SUBTYPE={C:'#13B8A6',A1:'#38BDF8',D:'#8B7BE8',G:'#173B4D',CRF02_AG:'#F59E0B',B:'#E5566D',F1:'#7ACBBF',CRF01_AE:'#B7ABF2'};
function treeViewer(){
  var host=$('#tv');if(!host)return;
  var treeOpts=(D.trees||[]).map(function(t,i){return '<option value="'+i+'">'+esc(t.label)+'</option>'}).join('');
  host.innerHTML='<h3 style="margin:.1rem 0 .3rem;font-size:1rem">Tree viewer</h3>'+
    '<p class="note" style="margin-top:0">Draws an Auspice tree here, with no download. Scroll to zoom the rows, Shift+scroll to zoom the branches, drag to pan, click a tip to pin its details. Internal branches are coloured only for reconstructed traits.</p>'+
    '<div class="tvbar">'+(treeOpts?'<div class="ctl"><label for="tvTree">Tree</label><select id="tvTree">'+treeOpts+'</select></div><div class="ctl"><button class="btn" id="tvLoad" type="button">Load tree</button></div>':'')+
    '<div class="ctl"><label for="tvFile">'+(treeOpts?'or a file':'Auspice JSON file')+'</label><input id="tvFile" type="file" accept=".json,application/json"></div>'+
    '<div class="ctl"><label for="tvColor">Colour by</label><select id="tvColor" disabled></select></div>'+
    '<div class="ctl"><span class="lab">Axis</span><div class="seg" id="tvAxis" role="group" aria-label="Axis"><button type="button" data-a="div" aria-pressed="true">Divergence</button><button type="button" data-a="date" aria-pressed="false">Date</button></div></div>'+
    '<div class="ctl"><label for="tvSearch">Find a tip</label><input id="tvSearch" type="search" placeholder="accession"></div>'+
    '<div class="ctl"><button class="btn ghost" id="tvFit" type="button">Fit</button></div></div>'+
    '<div class="tvwrap"><canvas id="tvCanvas" role="img" aria-label="Phylogenetic tree"></canvas><div id="tvTip"></div></div>'+
    '<div class="legend" id="tvLegend"></div><div id="tvPin"></div><p class="note" id="tvMsg">No tree loaded.</p>';
  var cv=$('#tvCanvas'),ctx=cv.getContext('2d'),tip=$('#tvTip');
  var V={nodes:[],leaves:[],meta:null,axis:'div',key:null,cs:null,hl:null,found:null,pin:null,kx:1,ky:1,tx:0,ty:0,xr:[0,1],hasDate:false,hover:null};
  var M={l:14,r:14,t:10,b:30};
  function msg(t){$('#tvMsg').textContent=t}
  function val(n,k){var a=n.node_attrs&&n.node_attrs[k];return a==null?undefined:(typeof a==='object'?a.value:a)}

  function setTree(json){
    var root=json.tree||json,nodes=[],leaves=[],st=[root];root._p=null;
    while(st.length){var n=st.pop();nodes.push(n);var ch=n.children;
      if(ch&&ch.length){for(var i=ch.length-1;i>=0;i--){ch[i]._p=n;st.push(ch[i])}}else{n._leaf=leaves.length;leaves.push(n)}}
    // rows: leaves in order, internal nodes at the mean of first and last child (children come before parents when reversed)
    for(var j=nodes.length-1;j>=0;j--){var m=nodes[j];
      if(m._leaf!=null)m._y=m._leaf;else m._y=(m.children[0]._y+m.children[m.children.length-1]._y)/2;
      m._xd=val(m,'div');var d=m.node_attrs&&m.node_attrs.num_date;m._xt=d==null?undefined:(typeof d==='object'?d.value:d)}
    nodes.forEach(function(m){var p=m._p;if(m._xd==null)m._xd=p?p._xd:0;if(m._xt==null)m._xt=p?p._xt:undefined});
    V.nodes=nodes;V.leaves=leaves;V.meta=json.meta||{};V.found=null;V.pin=null;V.hl=null;
    V.hasDate=leaves.some(function(l){return l._xt!=null});
    var segBtns=$('#tvAxis').querySelectorAll('button');segBtns[1].disabled=!V.hasDate;
    V.axis=V.hasDate&&/time tree/i.test(V.meta.title||'')?'date':'div';
    syncAxis();
    var cols=(V.meta.colorings||[]).filter(function(c){return leaves.some(function(l){return val(l,c.key)!=null})});
    var sel=$('#tvColor');sel.innerHTML=cols.map(function(c){return '<option value="'+esc(c.key)+'">'+esc(c.title||c.key)+'</option>'}).join('');sel.disabled=!cols.length;
    var def=(V.meta.display_defaults||{}).color_by;V.key=cols.some(function(c){return c.key===def})?def:(cols[0]&&cols[0].key);if(V.key)sel.value=V.key;
    buildScale();fit();
    msg((V.meta.title||'Tree')+': '+fmt(leaves.length)+' tips.');
  }
  function syncAxis(){$('#tvAxis').querySelectorAll('button').forEach(function(b){b.setAttribute('aria-pressed',String(b.dataset.a===V.axis))});
    var f=V.axis==='date'?'_xt':'_xd',lo=Infinity,hi=-Infinity;
    V.nodes.forEach(function(n){var x=n[f];if(x!=null){if(x<lo)lo=x;if(x>hi)hi=x}});
    if(!(hi>lo)){lo=0;hi=1}V.xr=[lo,hi];V.xf=f}
  function buildScale(){
    var k=V.key;if(!k){V.cs=null;renderLegend();return}
    var meta=(V.meta.colorings||[]).filter(function(c){return c.key===k})[0]||{},counts={};
    V.leaves.forEach(function(l){var v=val(l,k);if(v!=null){counts[v]=(counts[v]||0)+1}});
    var keys=Object.keys(counts);
    if(meta.type==='ordinal'||meta.type==='continuous'){
      var nums=keys.map(Number).filter(function(x){return !isNaN(x)}),lo=Math.min.apply(null,nums),hi=Math.max.apply(null,nums);
      V.cs={ord:true,lo:lo,hi:hi,color:function(v){var x=Number(v);if(isNaN(x))return '#CBD5E1';return lerp('#BFEFE9','#101D2D',hi>lo?(x-lo)/(hi-lo):.5)}};
    }else{
      keys.sort(function(a,b){return counts[b]-counts[a]});var map={},pi=0,fixedUsed=0;
      keys.forEach(function(v,i){
        if(k==='subtype'&&TV_SUBTYPE[v])map[v]=TV_SUBTYPE[v];
        else if(TV_FIXED[String(v).toLowerCase()])map[v]=TV_FIXED[String(v).toLowerCase()];
        else if(pi<TV_PALETTE.length){map[v]=TV_PALETTE[pi++]}
        else map[v]='#94A3B8'});
      V.cs={ord:false,counts:counts,keys:keys,map:map,other:keys.filter(function(v){return map[v]==='#94A3B8'}).length,
        color:function(v){return map[v]||'#CBD5E1'}};
    }
    renderLegend();
  }
  function renderLegend(){
    var el=$('#tvLegend'),cs=V.cs;if(!cs){el.innerHTML='';return}
    if(cs.ord){el.innerHTML='<span>'+esc(cs.lo)+'</span><span class="ramp" style="width:180px;background:linear-gradient(90deg,#BFEFE9,#101D2D)"></span><span>'+esc(cs.hi)+'</span>';return}
    var shown=cs.keys.filter(function(v){return cs.map[v]!=='#94A3B8'}).slice(0,16);
    el.innerHTML=shown.map(function(v){return '<button type="button" data-v="'+esc(v)+'" aria-pressed="'+(V.hl===v)+'"><i style="background:'+cs.map[v]+'"></i>'+esc(v)+' <span class="muted">'+fmt(cs.counts[v])+'</span></button>'}).join('')+(cs.other?'<span class="muted" style="font-size:.8rem">+ '+cs.other+' more values in grey</span>':'');
  }
  $('#tvLegend').addEventListener('click',function(e){var b=e.target.closest('button');if(!b)return;V.hl=V.hl===b.dataset.v?null:b.dataset.v;renderLegend();draw()});
  $('#tvColor').addEventListener('change',function(e){V.key=e.target.value;V.hl=null;buildScale();draw()});
  $('#tvAxis').addEventListener('click',function(e){var b=e.target.closest('button');if(!b||b.disabled)return;V.axis=b.dataset.a;syncAxis();draw()});
  $('#tvFit').addEventListener('click',function(){fit()});
  $('#tvSearch').addEventListener('input',function(e){
    var q=e.target.value.trim().toLowerCase();V.found=null;
    if(q.length>=3){var hit=V.leaves.filter(function(l){return String(l.name).toLowerCase().indexOf(q)>=0});
      V.found=hit.length?hit:null;msg(hit.length?fmt(hit.length)+' tips match.':'No tip matches.');
      if(hit.length===1){var H=cv.clientHeight,rowH=(H-M.t-M.b)/V.leaves.length;V.ky=Math.max(V.ky,Math.min(40,24/rowH));V.ty=(H/2)-(M.t+hit[0]._y*rowH*V.ky)}}
    draw()});
  function fit(){V.kx=1;V.ky=1;V.tx=0;V.ty=0;resize();}
  function resize(){var r=cv.getBoundingClientRect(),dpr=window.devicePixelRatio||1;cv.width=Math.round(r.width*dpr);cv.height=Math.round(r.height*dpr);draw()}
  window.addEventListener('resize',function(){if(V.nodes.length)resize()});

  function geom(){var W=cv.clientWidth,H=cv.clientHeight,rowH=(H-M.t-M.b)/Math.max(1,V.leaves.length),pw=W-M.l-M.r;
    return {W:W,H:H,rowH:rowH,pw:pw,
      sx:function(x){return M.l+(x-V.xr[0])/(V.xr[1]-V.xr[0])*pw*V.kx+V.tx},
      ix:function(px){return V.xr[0]+(px-V.tx-M.l)/(pw*V.kx)*(V.xr[1]-V.xr[0])},
      sy:function(y){return M.t+y*rowH*V.ky+V.ty}}}
  function niceStep(span,n){var raw=span/n,p=Math.pow(10,Math.floor(Math.log10(raw))),f=raw/p;return (f<1.5?1:f<3?2:f<7?5:10)*p}
  function draw(){
    if(!V.nodes.length){var c0=cv.getContext('2d'),dpr0=window.devicePixelRatio||1;c0.setTransform(dpr0,0,0,dpr0,0,0);c0.clearRect(0,0,cv.clientWidth,cv.clientHeight);c0.fillStyle='#64748B';c0.font='14px sans-serif';c0.fillText('Load a tree to see it here.',20,34);return}
    var g=geom(),dpr=cv.width/cv.clientWidth,f=V.xf;
    ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,g.W,g.H);ctx.lineCap='round';
    var rowPx=g.rowH*V.ky,cs=V.cs,dim=V.hl!=null,key=V.key;
    function colorOf(n){if(!cs||!key)return '#8FA1B3';var v=val(n,key);return v==null?null:cs.color(v)}
    function match(n){if(!dim)return true;return String(val(n,key))===V.hl}
    // pass 1: structure, neutral
    ctx.strokeStyle='#A7B5C3';ctx.lineWidth=1;ctx.beginPath();
    var nodes=V.nodes;
    for(var i=0;i<nodes.length;i++){var n=nodes[i];if(n._leaf!=null)continue;
      var x=g.sx(n[f]),y1=g.sy(n.children[0]._y),y2=g.sy(n.children[n.children.length-1]._y);
      if(y2<-4||y1>g.H+4)continue;
      ctx.moveTo(x,y1);ctx.lineTo(x,y2);
      for(var c=0;c<n.children.length;c++){var ch=n.children[c];if(ch._leaf!=null)continue;var cy=g.sy(ch._y);if(cy<-4||cy>g.H+4)continue;ctx.moveTo(x,cy);ctx.lineTo(g.sx(ch[f]),cy)}}
    ctx.stroke();
    // reconstructed-trait colours on internal branches
    if(key){for(var i2=0;i2<nodes.length;i2++){var q=nodes[i2];if(q._leaf!=null||!q._p)continue;var cq=colorOf(q);if(!cq)continue;var yq=g.sy(q._y);if(yq<-4||yq>g.H+4)continue;
      ctx.strokeStyle=cq;ctx.lineWidth=1.6;ctx.beginPath();ctx.moveTo(g.sx(q._p[f]),yq);ctx.lineTo(g.sx(q[f]),yq);ctx.stroke()}}
    // pass 2: tips
    var lw=Math.max(1,Math.min(3,rowPx*.8));
    for(var k=0;k<V.leaves.length;k++){var l=V.leaves[k],y=g.sy(l._y);if(y<-6||y>g.H+6)continue;
      var col=colorOf(l)||'#8FA1B3',on=match(l);ctx.globalAlpha=on?1:.12;
      var x1=g.sx(l._p?l._p[f]:l[f]),x2=g.sx(l[f]);
      ctx.strokeStyle=col;ctx.lineWidth=lw;ctx.beginPath();ctx.moveTo(x1,y);ctx.lineTo(x2,y);ctx.stroke();
      if(rowPx>=3){ctx.fillStyle=col;ctx.beginPath();ctx.arc(x2,y,Math.min(4.5,rowPx*.38),0,6.2832);ctx.fill()}}
    ctx.globalAlpha=1;
    if(V.found){V.found.slice(0,300).forEach(function(l){var y=g.sy(l._y),x=g.sx(l[f]);ctx.strokeStyle='#101D2D';ctx.lineWidth=2;ctx.beginPath();ctx.arc(x,y,7,0,6.2832);ctx.stroke()})}
    if(V.pin){var pl=V.pin,py=g.sy(pl._y),px=g.sx(pl[f]);ctx.strokeStyle='#101D2D';ctx.lineWidth=2.2;ctx.beginPath();ctx.arc(px,py,8,0,6.2832);ctx.stroke()}
    // axis
    ctx.fillStyle='#FFFFFF';ctx.fillRect(0,g.H-M.b+4,g.W,M.b);ctx.strokeStyle='#DCE5EC';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(0,g.H-M.b+4.5);ctx.lineTo(g.W,g.H-M.b+4.5);ctx.stroke();
    var lo=g.ix(M.l),hi=g.ix(g.W-M.r),step=niceStep(hi-lo,6),first=Math.ceil(lo/step)*step;
    ctx.fillStyle='#64748B';ctx.font='11px sans-serif';ctx.textAlign='center';
    for(var t=first;t<=hi+1e-9;t+=step){var tx=g.sx(t);ctx.beginPath();ctx.moveTo(tx,g.H-M.b+4);ctx.lineTo(tx,g.H-M.b+9);ctx.stroke();
      ctx.fillText(V.axis==='date'?String(Math.round(t*10)/10):String(Math.round(t*1e4)/1e4),tx,g.H-6)}
    ctx.textAlign='left';
  }
  // interaction
  var drag=null;
  cv.addEventListener('wheel',function(e){if(!V.nodes.length)return;e.preventDefault();
    var r=cv.getBoundingClientRect(),mx=e.clientX-r.left,my=e.clientY-r.top,z=Math.exp(-e.deltaY*0.0015);
    if(e.shiftKey||Math.abs(e.deltaX)>Math.abs(e.deltaY)){var nk=Math.max(1,Math.min(60,V.kx*z));V.tx=mx-(mx-V.tx)*(nk/V.kx);V.kx=nk}
    else{var nky=Math.max(1,Math.min(400,V.ky*z));V.ty=my-(my-V.ty)*(nky/V.ky);V.ky=nky}
    clamp();draw()},{passive:false});
  function clamp(){var g=geom();V.ty=Math.min(0,Math.max(V.ty,g.H-M.b-M.t-(g.H-M.t-M.b)*V.ky-20));V.tx=Math.min(0,Math.max(V.tx,-(g.pw*(V.kx-1))-20))}
  cv.addEventListener('pointerdown',function(e){if(!V.nodes.length)return;drag={x:e.clientX,y:e.clientY,tx:V.tx,ty:V.ty,moved:false};cv.setPointerCapture(e.pointerId);cv.classList.add('drag')});
  cv.addEventListener('pointermove',function(e){
    var r=cv.getBoundingClientRect(),mx=e.clientX-r.left,my=e.clientY-r.top;
    if(drag){var dx=e.clientX-drag.x,dy=e.clientY-drag.y;if(Math.abs(dx)+Math.abs(dy)>3)drag.moved=true;V.tx=drag.tx+dx;V.ty=drag.ty+dy;clamp();draw();tip.style.display='none';return}
    hoverAt(mx,my)});
  cv.addEventListener('pointerup',function(e){cv.classList.remove('drag');var was=drag;drag=null;
    if(was&&!was.moved){var r=cv.getBoundingClientRect(),l=tipAt(e.clientX-r.left,e.clientY-r.top);V.pin=l;pinPanel(l);draw()}});
  cv.addEventListener('pointerleave',function(){tip.style.display='none'});
  function tipAt(mx,my){if(!V.nodes.length)return null;var g=geom(),row=Math.round((my-M.t-V.ty)/(g.rowH*V.ky));
    if(row<0||row>=V.leaves.length)return null;var l=V.leaves[row],dy=Math.abs(g.sy(l._y)-my);
    if(dy>Math.max(7,g.rowH*V.ky/2+2))return null;if(mx<g.sx(l[V.xf])-55)return null;return l}
  function attrsOf(l){var o=[];['subtype','country','year','study','treatment_context','sdrm_mutations','cluster'].forEach(function(k){var v=val(l,k);if(v!=null)o.push([k.replace(/_/g,' '),v])});return o}
  function hoverAt(mx,my){var l=tipAt(mx,my);
    if(!l){tip.style.display='none';return}
    tip.innerHTML='<b>'+esc(l.name)+'</b><br>'+attrsOf(l).slice(0,5).map(function(a){return esc(a[0])+': '+esc(a[1])}).join('<br>');
    tip.style.display='block';var w=cv.clientWidth;tip.style.left=Math.min(w-310,mx+14)+'px';tip.style.top=Math.max(4,my-10)+'px'}
  function pinPanel(l){var el=$('#tvPin');if(!l){el.innerHTML='';return}
    var muts=l.branch_attrs&&l.branch_attrs.mutations,nn=muts&&muts.nuc?muts.nuc.length:null;
    el.innerHTML='<div class="callout"><strong>'+esc(l.name)+'</strong><table><tbody>'+attrsOf(l).map(function(a){return '<tr><td class="muted" style="width:11rem">'+esc(a[0])+'</td><td>'+esc(a[1])+'</td></tr>'}).join('')+
      (V.hasDate&&l._xt!=null?'<tr><td class="muted">date</td><td>'+esc(Math.round(l._xt*100)/100)+'</td></tr>':'')+(nn!=null?'<tr><td class="muted">nucleotide changes on the terminal branch</td><td>'+nn+'</td></tr>':'')+'</tbody></table></div>'}

  // loading
  function parse(text,label){try{setTree(JSON.parse(text))}catch(err){msg('Could not read '+label+' as an Auspice JSON.')}}
  var loadBtn=$('#tvLoad');
  if(loadBtn)loadBtn.addEventListener('click',function(){var t=D.trees[+$('#tvTree').value];msg('Loading '+t.name+'…');
    fetch(t.file).then(function(r){if(!r.ok)throw new Error('HTTP '+r.status);return r.text()}).then(function(x){parse(x,t.name)})
      .catch(function(){msg('Could not fetch '+t.file+'. Browsers block this when the page is opened from a file; serve the folder over http (for example python -m http.server in the site folder) or choose the file below.')})});
  $('#tvFile').addEventListener('change',function(e){var f=e.target.files[0];if(!f)return;msg('Reading '+f.name+'…');var r=new FileReader();r.onload=function(){parse(r.result,f.name)};r.readAsText(f)});
  ['dragover','drop'].forEach(function(ev){cv.addEventListener(ev,function(e){e.preventDefault();if(ev==='drop'&&e.dataTransfer.files[0]){var f=e.dataTransfer.files[0],r=new FileReader();r.onload=function(){parse(r.result,f.name)};r.readAsText(f)}})});
  resize();
}

/* ---------- limitations, references, footer ---------- */
$('#limList').innerHTML=D.limitations.map(function(t){return '<li>'+esc(t)+'</li>'}).join('');
$('#refList').innerHTML=D.references.map(function(t){return '<li>'+esc(t)+'</li>'}).join('');
var vs=Object.keys(D.list_versions).map(function(k){return (ROLE_LABEL[k]||k)+' '+D.list_versions[k]}).join('; ');
$('#foot').innerHTML='Virotrace explorer, scope <code>'+esc(D.scope)+'</code>, generated '+esc(D.generated)+'. Mutation lists: '+esc(vs)+'. Minimum group size '+esc(D.min_stratum_n)+'. Aggregate tables only: no accessions or per-sequence rows are in this page. The same data is in <code>data.json</code>.';

setRole(st.role);treeViewer();
})();
</script>
</body></html>
"""
