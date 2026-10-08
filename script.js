const $=id=>document.getElementById(id);
function num(id){return Number($(id).value)||0}
function slabDesign(){
 const L=num("sL"),h=num("sh"),fc=num("sfc"),fy=num("sfy"),cover=num("scover"),db=num("sbar");
 const finish=num("sfinish"),live=num("slive"),other=num("sother"),coef=num("support"),vcoef=num("vcoef");
 const sw=h/1000*24,D=sw+finish+other,wu=Math.max(1.4*D,1.2*D+1.6*live);
 const Mu=wu*L*L*coef,b=1000,d=h-cover-db/2,phi=.90;
 const R=Mu*1e6/(phi*b*d*d),disc=Math.max(0,1-2*R/(0.85*fc));
 let As=disc>0 ? (0.85*fc*b/fy)*(d-d*Math.sqrt(disc)) : 1e9;
 const Asmin=.0018*b*h,Asgov=Math.max(As,Asmin),area=Math.PI*db*db/4;
 const sReq=b*area/Asgov,sMax=Math.min(3*h,450);
 const s=Math.max(25,Math.floor(Math.min(sReq,sMax)/5)*5),Aprov=b*area/s;
 const Vu=wu*L*vcoef,Vc=2*Math.sqrt(fc)*b*d/1000,phiVc=.75*Vc;
 $("ssw").textContent=sw.toFixed(2);$("sd").textContent=D.toFixed(2);$("swu").textContent=wu.toFixed(2);
 $("smu").textContent=Mu.toFixed(2);$("sddep").textContent=d.toFixed(2);$("sasreq").textContent=As.toFixed(2);
 $("sasmin").textContent=Asmin.toFixed(2);$("sasgov").textContent=Asgov.toFixed(2);
 $("sprovide").textContent="Ø"+db+" @ "+s+" mm c/c ("+Aprov.toFixed(0)+" mm²/m)";
 $("svu").textContent=Vu.toFixed(2);$("svc").textContent=phiVc.toFixed(2);
 $("svstatus").textContent=phiVc>=Vu?"OK — shear capacity adequate":"NOT OK — increase thickness / redesign";
}

function solveLinear(A,b){
 const n=b.length,M=A.map((r,i)=>r.slice().concat([b[i]]));
 for(let k=0;k<n;k++){
  let p=k;for(let i=k+1;i<n;i++)if(Math.abs(M[i][k])>Math.abs(M[p][k]))p=i;
  if(Math.abs(M[p][k])<1e-12)throw new Error("Singular stiffness matrix — check supports.");
  [M[k],M[p]]=[M[p],M[k]];
  for(let i=k+1;i<n;i++){const q=M[i][k]/M[k][k];for(let j=k;j<=n;j++)M[i][j]-=q*M[k][j];}
 }
 const x=Array(n).fill(0);
 for(let i=n-1;i>=0;i--){let s=M[i][n];for(let j=i+1;j<n;j++)s-=M[i][j]*x[j];x[i]=s/M[i][i];}
 return x;
}
function fmtMatrix(A,d=3){return A.map(r=>r.map(v=>Number(v).toFixed(d)).join("   ")).join("\n")}

let lastBeamResult=null;
function beamDesign(){
 try{
  if(!lastBeamResult)beamStiffness();
  const fc=num("bfc"),fy=num("bfy"),cover=num("bcover"),db=num("bbar"),ds=num("bstir"),b=num("bb"),h=num("bh");
  const Mu=lastBeamResult.Mu,Vu=lastBeamResult.Vu,d=h-cover-ds-db/2,phi=.90;
  const R=Mu*1e6/(phi*b*d*d),disc=1-2*R/(0.85*fc);
  if(d<=0||disc<0)throw new Error("Section depth/material inputs are outside this preliminary design model.");
  const As=(0.85*fc*b/fy)*(d-d*Math.sqrt(disc));
  const Asmin=Math.max(1.4/fy*b*d,0.25*Math.sqrt(fc)/fy*b*d);
  const Asg=Math.max(As,Asmin);
  const area=Math.PI*db*db/4;
  let n=Math.max(2,Math.ceil(Asg/area)),Aprov=n*area;
  const Vc=0.17*Math.sqrt(fc)*b*d/1000,phiVc=.75*Vc;
  const Vs=Math.max(0,Vu/0.75-Vc);
  const Asv=2*Math.PI*ds*ds/4;
  let s=Vs>0 ? Asv*fy*d/(Vs*1000) : 9999;
  const smax=Math.min(d/2,600);
  s=Math.max(75,Math.floor(Math.min(s,smax)/25)*25);
  $("dMu").textContent=Mu.toFixed(2);$("dVu").textContent=Vu.toFixed(2);$("dd").textContent=d.toFixed(1);
  $("dAs").textContent=As.toFixed(0);$("dAsmin").textContent=Asmin.toFixed(0);
  $("dBars").textContent=n+"Ø"+db; $("dAsprov").textContent=Aprov.toFixed(0);
  $("dVc").textContent=phiVc.toFixed(2);
  $("dStir").textContent=Vs<=0?"Ø"+ds+" 2-leg @ "+smax.toFixed(0)+" mm max":"Ø"+ds+" 2-leg @ "+s+" mm c/c";
  $("dSmax").textContent=smax.toFixed(0);
 }catch(e){alert(e.message)}
}

function beamStiffness(){
 try{
  const L=num("bL"),b=num("bb")/1000,h=num("bh")/1000,E=num("bE")*1000,w=num("bw"),P=num("bP"),a=Math.min(L,Math.max(0,num("ba")));
  if(L<=0||b<=0||h<=0||E<=0)throw new Error("Check geometry and E.");
  const I=b*Math.pow(h,3)/12,EI=E*I;
  const c=EI/Math.pow(L,3),K=[
   [12*c,6*L*c,-12*c,6*L*c],
   [6*L*c,4*L*L*c,-6*L*c,2*L*L*c],
   [-12*c,-6*L*c,12*c,-6*L*c],
   [6*L*c,2*L*L*c,-6*L*c,4*L*L*c]
  ];
  const F=[-w*L/2,-w*L*L/12,-w*L/2,w*L*L/12];
  if(P>0){const xi=a/L,Fp=[-P*(1-3*xi*xi+2*xi*xi*xi),-P*L*xi*(1-xi)*(1-xi),-P*(3*xi*xi-2*xi*xi*xi),P*L*xi*xi*(1-xi)];for(let i=0;i<4;i++)F[i]+=Fp[i]}
  const fixed=[];
  if(numSupport("bLeft")==="fixed"){fixed.push(0,1)}else fixed.push(0);
  if(numSupport("bRight")==="fixed"){fixed.push(2,3)}else if(numSupport("bRight")==="roller"){fixed.push(2)}
  const free=[0,1,2,3].filter(i=>!fixed.includes(i)),Kr=free.map(i=>free.map(j=>K[i][j])),Fr=free.map(i=>F[i]);
  const d=Array(4).fill(0),dr=solveLinear(Kr,Fr);free.forEach((i,n)=>d[i]=dr[n]);
  const reactions=K.map((r,i)=>r.reduce((s,v,j)=>s+v*d[j],0)-F[i]);
  const r1=reactions[0],r2=reactions[2];
  const samples=80,pts=[];let maxM=0,maxV=0,maxDef=0,midDef=0;
  for(let i=0;i<=samples;i++){const x=L*i/samples,V=r1-w*x-(x>=a?P:0),M=r1*x-w*x*x/2-(x>=a?P*(x-a):0);pts.push({x,V,M});maxM=Math.max(maxM,Math.abs(M));maxV=Math.max(maxV,Math.abs(V));if(Math.abs(x-L/2)<L/samples)midDef=Math.abs(d[0]+(d[2]-d[0])*.5)}
  maxDef=Math.max(...pts.map((_,i)=>Math.abs(d[0]+(d[2]-d[0])*i/samples)));
  lastBeamResult={Mu:maxM,Vu:maxV};
  $("bvmax").textContent=(maxDef*1000).toFixed(3);
  $("br1").textContent=r1.toFixed(2);$("br2").textContent=r2.toFixed(2);
  $("bm1").textContent=pts[0].M.toFixed(2);$("bm2").textContent=pts[pts.length-1].M.toFixed(2);
  $("bMmax").textContent=maxM.toFixed(2);$("bVmax").textContent=maxV.toFixed(2);
  $("bK").textContent=fmtMatrix(K,2);$("bPvec").textContent=F.map(v=>v.toFixed(3)).join("\n");$("bD").textContent=d.map(v=>v.toFixed(8)).join("\n");
  let old=document.getElementById("bDiagram");if(old)old.remove();
  const wrap=document.createElement("div");wrap.className="diagram-card";wrap.id="bDiagram";
  wrap.innerHTML="<h2>Shear Force & Bending Moment</h2><p>Sign convention: positive shear upward on the left cut; positive sagging moment.</p><div class='charts'><div><h3>SFD</h3><canvas id='sfdCanvas' width='760' height='240'></canvas></div><div><h3>BMD</h3><canvas id='bmdCanvas' width='760' height='240'></canvas></div></div>";
  document.querySelector(".matrix-card").after(wrap);drawBeamChart("sfdCanvas",pts,"V","Shear (kN)");drawBeamChart("bmdCanvas",pts,"M","Moment (kN·m)");
 }catch(e){alert(e.message)}
}
function numSupport(id){return $(id).value}
function drawBeamChart(id,pts,key,label){
 const c=$(id),ctx=c.getContext("2d"),W=c.width,H=c.height,pad=42;
 ctx.clearRect(0,0,W,H);ctx.strokeStyle="#cfd8e4";ctx.lineWidth=1;
 const vals=pts.map(p=>p[key]),max=Math.max(...vals.map(Math.abs),1),x0=pad,x1=W-pad,y0=H/2,scale=(H*.38)/max;
 ctx.beginPath();ctx.moveTo(x0,y0);ctx.lineTo(x1,y0);ctx.stroke();
 ctx.strokeStyle="#1674e8";ctx.lineWidth=2;ctx.beginPath();
 pts.forEach((p,i)=>{const x=x0+(x1-x0)*i/(pts.length-1),y=y0-p[key]*scale;if(i)ctx.lineTo(x,y);else ctx.moveTo(x,y)});ctx.stroke();
 ctx.fillStyle="#66758a";ctx.font="11px system-ui";ctx.fillText("0",8,y0+4);ctx.fillText(max.toFixed(1),8,y0-scale+4);ctx.fillText((-max).toFixed(1),2,y0+scale+4);ctx.fillText(label,x0,18);ctx.fillText("0",x0-4,H-12);ctx.fillText(pts[pts.length-1].x.toFixed(2)+" m",x1-35,H-12);
}

function calculate(){
 const t=num("slab_t")/1000,gc=num("conc_w"),finish=num("finish"),live=num("occupancy");
 const bb=num("beam_b")/1000,bh=num("beam_h")/1000;
 const wt=num("wall_t")/1000,wh=num("wall_h"),ww=num("wall_w"),pl=num("plaster");
 const slab=t*gc,dead=slab+finish,service=dead+live;
 const beam=bb*bh*gc;
 const wall=wt*wh*ww+2*wh*pl;
 const u14=1.4*dead, u12=1.2*dead+1.6*live, factored=Math.max(u14,u12);
 $("slab_out").textContent=slab.toFixed(2);
 $("finish_out").textContent=finish.toFixed(2);
 $("dead_out").textContent=dead.toFixed(2);
 $("live_out").textContent=live.toFixed(2);
 $("service_out").textContent=service.toFixed(2);
 $("beam_out").textContent=beam.toFixed(2);
 $("wall_out").textContent=wall.toFixed(2);
 $("u14_out").textContent=u14.toFixed(2);
 $("u12_out").textContent=u12.toFixed(2);
 $("factored_out").textContent=factored.toFixed(2);
 $("live_selected").textContent=live.toFixed(2)+" kN/m²";
}

function drawBeamModel(){const c=$("beamModelCanvas");if(!c)return;const ctx=c.getContext("2d"),W=c.width,H=c.height;ctx.clearRect(0,0,W,H);ctx.fillStyle="#f7f9fc";ctx.fillRect(0,0,W,H);const L=Math.max(.1,num("bL")),x1=110,x2=W-110,y=270,s=(x2-x1)/L;ctx.strokeStyle="#cfd7e2";ctx.lineWidth=1;for(let x=x1;x<=x2;x+=s){ctx.beginPath();ctx.moveTo(x,70);ctx.lineTo(x,470);ctx.stroke()}ctx.strokeStyle="#27364b";ctx.lineWidth=8;ctx.beginPath();ctx.moveTo(x1,y);ctx.lineTo(x2,y);ctx.stroke();ctx.lineWidth=2;ctx.strokeStyle="#1674e8";ctx.beginPath();ctx.moveTo(x1,y-4);ctx.lineTo(x2,y-4);ctx.stroke();function support(x,type){ctx.strokeStyle="#334155";ctx.lineWidth=2;if(type==="fixed"){ctx.beginPath();ctx.moveTo(x,y-32);ctx.lineTo(x,y+32);ctx.stroke();for(let yy=y-28;yy<=y+28;yy+=9){ctx.beginPath();ctx.moveTo(x,yy);ctx.lineTo(x-13,yy+7);ctx.stroke()}}else{ctx.beginPath();ctx.moveTo(x,y+6);ctx.lineTo(x-18,y+34);ctx.lineTo(x+18,y+34);ctx.closePath();ctx.stroke();ctx.beginPath();ctx.moveTo(x-25,y+40);ctx.lineTo(x+25,y+40);ctx.stroke()}}support(x1,$("bLeft").value);support(x2,$("bRight").value);const w=num("bw");if(w>0){ctx.strokeStyle="#d14b4b";ctx.fillStyle="#d14b4b";ctx.lineWidth=1.5;for(let i=0;i<=10;i++){const x=x1+(x2-x1)*i/10;ctx.beginPath();ctx.moveTo(x,y-55);ctx.lineTo(x,y-12);ctx.stroke();ctx.beginPath();ctx.moveTo(x-5,y-20);ctx.lineTo(x,y-12);ctx.lineTo(x+5,y-20);ctx.fill()}}const P=num("bP"),a=Math.min(L,Math.max(0,num("ba")));if(P>0){const x=x1+a*s;ctx.strokeStyle="#7c3aed";ctx.fillStyle="#7c3aed";ctx.lineWidth=3;ctx.beginPath();ctx.moveTo(x,y-90);ctx.lineTo(x,y-18);ctx.stroke();ctx.beginPath();ctx.moveTo(x-7,y-28);ctx.lineTo(x,y-16);ctx.lineTo(x+7,y-28);ctx.fill()}ctx.fillStyle="#34435a";ctx.font="12px system-ui";ctx.fillText("Beam 1",x1,y-105);ctx.fillText("L = "+L.toFixed(2)+" m",(x1+x2)/2-30,y+70);ctx.fillText($("bLeft").value.toUpperCase(),x1-25,y+62);ctx.fillText($("bRight").value.toUpperCase(),x2-25,y+62);}
function updateCad(){drawBeamModel();$("cadB").textContent=num("bb").toFixed(0)+" mm";$("cadH").textContent=num("bh").toFixed(0)+" mm";$("cadCoord").textContent="L = "+num("bL").toFixed(2)+" m";if(lastBeamResult){$("cadMu").textContent=lastBeamResult.Mu.toFixed(2)+" kN·m";$("cadVu").textContent=lastBeamResult.Vu.toFixed(2)+" kN";$("cadDef").textContent=$("bvmax").textContent+" mm"}}
document.querySelectorAll("#beam input,#beam select").forEach(i=>i.addEventListener("input",updateCad));
$("cadAnalyze").addEventListener("click",()=>{beamStiffness();updateCad()});$("cadDesign").addEventListener("click",()=>{beamDesign();document.getElementById("beamDesign").scrollIntoView({behavior:"smooth",block:"center"})});
setTimeout(updateCad,50);
$("calculate").addEventListener("click",calculate);
$("beamCalculate").addEventListener("click",beamStiffness);
function parseSeries(id,n,def=0){const raw=$(id).value.split(',').map(s=>Number(s.trim())).filter(v=>Number.isFinite(v));if(!raw.length)return Array(n).fill(def);const out=[];for(let i=0;i<n;i++)out.push(raw[Math.min(i,raw.length-1)]);return out}
function elementData(L,b,h,E,w,P,a){const I=b*Math.pow(h,3)/12,EI=E*I,c=EI/Math.pow(L,3),k=[[12*c,6*L*c,-12*c,6*L*c],[6*L*c,4*L*L*c,-6*L*c,2*L*L*c],[-12*c,-6*L*c,12*c,-6*L*c],[6*L*c,2*L*L*c,-6*L*c,4*L*L*c]],f=[-w*L/2,-w*L*L/12,-w*L/2,w*L*L/12];if(P){const x=Math.max(0,Math.min(L,a)),xi=x/L,fp=[-P*(1-3*xi*xi+2*xi*xi*xi),-P*L*xi*(1-xi)*(1-xi),-P*(3*xi*xi-2*xi*xi*xi),P*L*xi*xi*(1-xi)];for(let i=0;i<4;i++)f[i]+=fp[i]}return{k,f}}
function continuousBeamStiffness(){try{const n=Number($("cbN").value),b=num("cbB")/1000,h=num("cbH")/1000,E=num("cbE")*1000;if(n<2||n>4||b<=0||h<=0||E<=0)throw new Error("Check continuous-beam inputs.");const Ls=parseSeries("cbLengths",n,6).map(x=>Math.max(.2,x)),ws=parseSeries("cbLoads",n,10).map(x=>Math.max(0,x)),Ps=parseSeries("cbPoints",n,0).map(x=>Math.max(0,x)),As=parseSeries("cbPositions",n,0),N=n+1,dof=2*N,K=Array.from({length:dof},()=>Array(dof).fill(0)),F=Array(dof).fill(0),els=[];for(let e=0;e<n;e++){const L=Ls[e],a=Math.max(0,Math.min(L,As[e])),ed=elementData(L,b,h,E,ws[e],Ps[e],a),map=[2*e,2*e+1,2*e+2,2*e+3];els.push({...ed,L,w:ws[e],P:Ps[e],a});for(let i=0;i<4;i++){F[map[i]]+=ed.f[i];for(let j=0;j<4;j++)K[map[i]][map[j]]+=ed.k[i][j]}}const fixed=[0,2*N-2];if($("cbLeft").value==="fixed")fixed.push(1);if($("cbRight").value==="fixed")fixed.push(2*N-1);for(let node=1;node<n;node++)fixed.push(2*node);const free=Array.from({length:dof},(_,i)=>i).filter(i=>!fixed.includes(i)),Kr=free.map(i=>free.map(j=>K[i][j])),Fr=free.map(i=>F[i]),D=Array(dof).fill(0),Dr=solveLinear(Kr,Fr);free.forEach((i,k)=>D[i]=Dr[k]);const reactions=K.map((r,i)=>r.reduce((s,v,j)=>s+v*D[j],0)-F[i]),spanResults=[],globalPts=[];let xGlobal=0,maxPos=0,maxNeg=0,maxV=0,maxDef=0;for(let e=0;e<n;e++){const el=els[e],map=[2*e,2*e+1,2*e+2,2*e+3],de=map.map(i=>D[i]),q=el.k.map((r,i)=>r.reduce((s,v,j)=>s+v*de[j],0)-el.f[i]),pts=[],samples=60;for(let i=0;i<=samples;i++){const x=el.L*i/samples;let V=q[0]-el.w*x,M=q[1]+q[0]*x-el.w*x*x/2;if(el.P&&x>=el.a){V-=el.P;M-=el.P*(x-el.a)}pts.push({x:xGlobal+x,V,M});maxV=Math.max(maxV,Math.abs(V));if(M>=0)maxPos=Math.max(maxPos,M);else maxNeg=Math.min(maxNeg,M);globalPts.push({x:xGlobal+x,V,M})}spanResults.push({q,positive:Math.max(...pts.map(p=>p.M)),negative:Math.min(...pts.map(p=>p.M)),endLeft:q[1],endRight:-q[3]});xGlobal+=el.L}for(let i=0;i<N;i++)maxDef=Math.max(maxDef,Math.abs(D[2*i]));$("cbPosM").textContent=maxPos.toFixed(2);$("cbNegM").textContent=Math.abs(maxNeg).toFixed(2);$("cbMaxV").textContent=maxV.toFixed(2);$("cbMaxD").textContent=(maxDef*1000).toFixed(3);$("cbK").textContent=fmtMatrix(K,2);$("cbPvec").textContent=F.map(v=>v.toFixed(3)).join("\\n");$("cbD").textContent=D.map(v=>v.toFixed(8)).join("\\n");let table="<table class='action-table'><thead><tr><th>Span / Joint</th><th>Length</th><th>Left M</th><th>Right M</th><th>+M</th><th>−M</th></tr></thead><tbody>";spanResults.forEach((r,i)=>{table+="<tr><td>Span "+(i+1)+"</td><td>"+Ls[i].toFixed(2)+" m</td><td>"+r.endLeft.toFixed(2)+"</td><td>"+r.endRight.toFixed(2)+"</td><td>"+r.positive.toFixed(2)+"</td><td>"+Math.abs(r.negative).toFixed(2)+"</td></tr>"});table+="</tbody></table><div class='joint-grid'>";for(let i=0;i<N;i++){table+="<div class='joint-box'><b>Joint "+(i+1)+"</b><span>v = "+(D[2*i]*1000).toFixed(3)+" mm</span><span>θ = "+D[2*i+1].toExponential(3)+" rad</span><span>Vertical reaction = "+reactions[2*i].toFixed(2)+" kN</span></div>"}table+="</div>";$("cbTableWrap").innerHTML=table;drawBeamChart("cbSfdCanvas",globalPts,"V","Shear (kN)");drawBeamChart("cbBmdCanvas",globalPts,"M","Moment (kN·m)")}catch(e){alert(e.message)}}$("cbCalculate").addEventListener("click",continuousBeamStiffness);continuousBeamStiffness();$("beamDesign").addEventListener("click",beamDesign);
document.querySelectorAll("#beam input,#beam select").forEach(i=>i.addEventListener("input",beamStiffness));
$("slabCalculate").addEventListener("click",slabDesign);
document.querySelectorAll("#slab input,#slab select").forEach(i=>i.addEventListener("input",slabDesign));
document.querySelectorAll("input,select").forEach(i=>i.addEventListener("input",calculate));
document.querySelectorAll(".nav-item").forEach(btn=>btn.addEventListener("click",()=>{
 document.querySelectorAll(".nav-item").forEach(x=>x.classList.remove("active"));btn.classList.add("active");
 document.querySelectorAll(".section,.placeholder").forEach(x=>x.style.display="none");
 const t=$(btn.dataset.section);if(t)t.style.display="block";
}));
calculate();
/* Graphical model builder — current analysis engine supports the active two-joint member */
let cadNodes=[
 {id:1,x:0,support:"pin"},
 {id:2,x:null,support:"roller"}
];
let cadPointLoads=[];
let cadUdl=10;

function syncGraphModel(){
 const L=Math.max(.5,num("bL"));
 if(cadNodes.length<2)cadNodes=[{id:1,x:0,support:"pin"},{id:2,x:L,support:"roller"}];
 if(cadNodes[1].x===null||cadNodes[1].x<=0)cadNodes[1].x=L;
 $("bL").value=cadNodes[1].x.toFixed(2);
 $("bLeft").value=cadNodes[0].support==="fixed"?"fixed":"pin";
 $("bRight").value=cadNodes[1].support==="fixed"?"fixed":cadNodes[1].support==="free"?"free":"roller";
 $("bw").value=cadUdl;
 if(cadPointLoads.length){
   const p=cadPointLoads[cadPointLoads.length-1];
   $("bP").value=p.P;$("ba").value=p.x.toFixed(2);
 }
}
function graphCanvasGeometry(){
 const c=$("beamModelCanvas"),L=Math.max(.5,num("bL")),x1=110,x2=c.width-110,y=270;
 return {c,L,x1,x2,y,s:(x2-x1)/L};
}
function graphHitNode(p){
 for(const n of cadNodes){
   const nx=p.x1+n.x*p.s;
   if(Math.abs(p.px-nx)<28 && Math.abs(p.py-p.y)<45)return n;
 }
 return null;
}
function drawGraphModel(){
 const {c,L,x1,x2,y,s}=graphCanvasGeometry(),ctx=c.getContext("2d"),W=c.width,H=c.height;
 ctx.clearRect(0,0,W,H);ctx.fillStyle="#f7f9fc";ctx.fillRect(0,0,W,H);
 ctx.strokeStyle="#d7dee8";ctx.lineWidth=1;
 for(let x=x1;x<=x2+1;x+=Math.max(45,s)){ctx.beginPath();ctx.moveTo(x,70);ctx.lineTo(x,470);ctx.stroke()}
 const pts=cadNodes.slice().sort((a,b)=>a.x-b.x);
 ctx.strokeStyle="#27364b";ctx.lineWidth=8;
 for(let i=0;i<pts.length-1;i++){const xa=x1+pts[i].x*s,xb=x1+pts[i+1].x*s;ctx.beginPath();ctx.moveTo(xa,y);ctx.lineTo(xb,y);ctx.stroke()}
 function support(x,type){
   ctx.strokeStyle="#334155";ctx.lineWidth=2;
   if(type==="fixed"){ctx.beginPath();ctx.moveTo(x,y-34);ctx.lineTo(x,y+34);ctx.stroke();
     for(let yy=y-30;yy<=y+30;yy+=9){ctx.beginPath();ctx.moveTo(x,yy);ctx.lineTo(x-14,yy+7);ctx.stroke()}
   }else if(type!=="free"){ctx.beginPath();ctx.moveTo(x,y+7);ctx.lineTo(x-18,y+34);ctx.lineTo(x+18,y+34);ctx.closePath();ctx.stroke();
     ctx.beginPath();ctx.moveTo(x-25,y+40);ctx.lineTo(x+25,y+40);ctx.stroke()}
 }
 for(const n of pts){const x=x1+n.x*s;support(x,n.support);ctx.fillStyle=n.id===cadSelection?"#1674e8":"#263548";ctx.beginPath();ctx.arc(x,y,7,0,Math.PI*2);ctx.fill();
   ctx.fillStyle="#34435a";ctx.font="11px system-ui";ctx.fillText("J"+n.id,x-8,y-48);ctx.fillText(n.support.toUpperCase(),x-25,y+62)}
 if(cadUdl>0&&pts.length>=2){ctx.strokeStyle="#d14b4b";ctx.fillStyle="#d14b4b";ctx.lineWidth=1.5;
   const xa=x1+pts[0].x*s,xb=x1+pts[pts.length-1].x*s;
   for(let i=0;i<=12;i++){const x=xa+(xb-xa)*i/12;ctx.beginPath();ctx.moveTo(x,y-55);ctx.lineTo(x,y-12);ctx.stroke();ctx.beginPath();ctx.moveTo(x-5,y-20);ctx.lineTo(x,y-12);ctx.lineTo(x+5,y-20);ctx.fill()}
   ctx.fillText("UDL = "+cadUdl.toFixed(2)+" kN/m",xa,y-68)
 }
 for(const p of cadPointLoads){const x=x1+p.x*s;ctx.strokeStyle="#7c3aed";ctx.fillStyle="#7c3aed";ctx.lineWidth=3;ctx.beginPath();ctx.moveTo(x,y-100);ctx.lineTo(x,y-18);ctx.stroke();ctx.beginPath();ctx.moveTo(x-7,y-28);ctx.lineTo(x,y-16);ctx.lineTo(x+7,y-28);ctx.fill();ctx.font="11px system-ui";ctx.fillText(p.P+" kN",x+8,y-80)}
 ctx.fillStyle="#1674e8";ctx.font="bold 11px system-ui";ctx.fillText("GRAPHICAL MODEL • "+cadMode.toUpperCase(),x1,32);
 ctx.fillStyle="#66758a";ctx.font="11px system-ui";ctx.fillText("Click empty space = create Joint",x1,500);
 ctx.fillText("Click Joint = select / support",x1+210,500);
 ctx.fillText("Loads mode: click beam = Point Load",x1+440,500);
}
function rebuildGraphTree(){
 const ex=document.querySelector(".model-explorer");if(!ex)return;
 ex.innerHTML='<div class="explorer-title">MODEL EXPLORER • '+cadMode.toUpperCase()+'</div><div class="tree-root">▾ <b>Beam Model</b></div><div class="tree-item selected">▾ Geometry</div>';
 cadNodes.forEach(n=>ex.insertAdjacentHTML("beforeend",'<div class="tree-sub cad-tree-node" data-node="'+n.id+'">Joint '+n.id+' — '+n.support+'</div>'));
 ex.insertAdjacentHTML("beforeend",'<div class="tree-sub cad-tree-node">Beam 1</div><div class="tree-item">▾ Loads</div><div class="tree-sub">UDL — '+cadUdl.toFixed(2)+' kN/m</div>');
 cadPointLoads.forEach((p,i)=>ex.insertAdjacentHTML("beforeend",'<div class="tree-sub">Point Load '+(i+1)+' — '+p.P+' kN @ '+p.x.toFixed(2)+' m</div>'));
 ex.insertAdjacentHTML("beforeend",'<div class="tree-item">▾ Analysis</div><div class="tree-sub">Direct Stiffness</div><div class="tree-item">▾ Design</div><div class="tree-sub">Flexure</div><div class="tree-sub">Shear</div><div class="tree-sub">Detailing</div>');
 ex.querySelectorAll(".cad-tree-node").forEach(el=>el.addEventListener("click",()=>{
   const n=cadNodes.find(q=>"Joint "+q.id+" — "+q.support===el.textContent.trim());
   if(n){cadSelection=n.id;showGraphProperties(n);drawGraphModel()}
 }));
}
function showGraphProperties(n){
 const panel=document.querySelector(".property-panel");if(!panel)return;
 const rows=panel.querySelectorAll(".prop-row");if(rows[0])rows[0].querySelector("strong").textContent="Joint";
 if(rows[1])rows[1].querySelector("strong").textContent="Joint "+n.id;
 if(rows[2])rows[2].querySelector("strong").textContent=(n.x||0).toFixed(2)+" m";
}
function addGraphJoint(x){
 const L=Math.max(.5,num("bL"));
 if(cadNodes.length>=4){alert("This graphical builder currently supports up to 4 joints.");return}
 const nx=Math.max(.25,Math.min(L-.25,x));
 if(cadNodes.some(n=>Math.abs(n.x-nx)<.25))return;
 const id=Math.max(...cadNodes.map(n=>n.id))+1;
 cadNodes.push({id,x:nx,support:"free"});
 cadNodes.sort((a,b)=>a.x-b.x);
 rebuildGraphTree();drawGraphModel();
}
function cycleSupport(n){
 const order=["free","pin","roller","fixed"],i=order.indexOf(n.support);n.support=order[(i+1)%order.length];
 if(n.id===1)$("bLeft").value=n.support==="fixed"?"fixed":"pin";
 if(n.id===cadNodes[cadNodes.length-1].id)$("bRight").value=n.support==="fixed"?"fixed":"roller";
 rebuildGraphTree();drawGraphModel();
}
function graphClick(ev){
 const p=cadCanvasPoint(ev),node=graphHitNode(p);
 if(cadMode==="model"){
   if(node){cadSelection=node.id;cycleSupport(node);showGraphProperties(node)}
   else addGraphJoint(p.x);
   syncGraphModel();drawGraphModel();
 }else if(cadMode==="loads"){
   const gg=graphCanvasGeometry();
   if(p.py>gg.y-100 && p.py<gg.y+45 && p.x>=0&&p.x<=gg.L){
     const P=Number(prompt("Point load P (kN)",num("bP")||10));
     if(Number.isFinite(P)&&P>0){cadPointLoads.push({P,x:p.x});$("bP").value=P;$("ba").value=p.x.toFixed(2);beamStiffness()}
     rebuildGraphTree();drawGraphModel();
   }
 }else if(cadMode==="analyze"){beamStiffness();setCadMode("results")}
 else if(cadMode==="design"){beamDesign();drawGraphModel()}
}
drawBeamModel=function(){syncGraphModel();drawGraphModel();};
const cadCanvas=$("beamModelCanvas");
if(cadCanvas){cadCanvas.addEventListener("click",graphClick)}
document.querySelectorAll(".cad-toolbar .tool[data-cad]").forEach(b=>b.addEventListener("click",()=>{setCadMode(b.dataset.cad);rebuildGraphTree()}));
const modelBtn=document.querySelector(".cad-toolbar .tool[data-cad='model']");
if(modelBtn)modelBtn.title="Model: click empty canvas to create joints; click joint to cycle support";
rebuildGraphTree();drawGraphModel();

/* Graphical model -> multi-span direct stiffness bridge */
function analyzeGraphicalModel(){
 try{
  cadNodes.sort((a,b)=>a.x-b.x);
  if(cadNodes.length<2)throw new Error("Create at least two joints.");
  const spans=cadNodes.length-1,Ls=cadNodes.slice(1).map((n,i)=>n.x-cadNodes[i].x);
  if(Ls.some(x=>x<=0.1))throw new Error("Joint spacing must be greater than 0.10 m.");
  const E=num("cbE")||num("bE")||25000,b=num("cbB")||num("bb"),h=num("cbH")||num("bh");
  const w=cadUdl;
  const pointBySpan=Array.from({length:spans},()=>[]);
  cadPointLoads.forEach(p=>{
   let k=Ls.findIndex((L,i)=>p.x>=cadNodes[i].x && p.x<=cadNodes[i+1].x);
   if(k<0)k=spans-1;
   pointBySpan[k].push({P:p.P,a:Math.max(0,Math.min(Ls[k],p.x-cadNodes[k].x))});
  });
  $("cbN").value=spans;$("cbB").value=b;$("cbH").value=h;$("cbE").value=E;
  $("cbLengths").value=Ls.map(x=>x.toFixed(2)).join(",");
  $("cbLoads").value=Ls.map(()=>w.toFixed(2)).join(",");
  $("cbPoints").value=Ls.map((_,i)=>pointBySpan[i][0]?.P||0).join(",");
  $("cbPositions").value=Ls.map((L,i)=>pointBySpan[i][0]?.a||0).join(",");
  const left=cadNodes[0].support,right=cadNodes[cadNodes.length-1].support;
  $("cbLeft").value=left==="fixed"?"fixed":"pin";
  $("cbRight").value=right==="fixed"?"fixed":"roller";
  continuousBeamStiffness();
  setCadMode("results");
  showGraphAnalysisSummary(Ls);
 }catch(e){alert(e.message)}
}
function showGraphAnalysisSummary(Ls){
 const panel=document.querySelector(".property-panel");if(!panel)return;
 let box=document.getElementById("graphAnalysisSummary");
 if(!box){box=document.createElement("div");box.id="graphAnalysisSummary";box.className="graph-summary";panel.appendChild(box)}
 box.innerHTML="<b>GRAPHICAL ANALYSIS</b><span>"+Ls.length+" span(s)</span><span>Total length: "+Ls.reduce((a,b)=>a+b,0).toFixed(2)+" m</span><span>UDL: "+cadUdl.toFixed(2)+" kN/m</span><span>Point loads: "+cadPointLoads.length+"</span><span>Model: Direct Stiffness</span>";
}
function graphicalResultsOverlay(){
 const c=$("beamModelCanvas");if(!c)return;
 const ctx=c.getContext("2d"),{L,x1,x2,y,s}=graphCanvasGeometry();
 if(cadMode!=="results"||!lastBeamResult)return;
 ctx.save();ctx.fillStyle="#1674e8";ctx.font="bold 11px system-ui";
 ctx.fillText("ANALYSIS RESULTS • SFD / BMD available below",x1,y+92);
 ctx.restore();
}
const oldGraphDraw=drawGraphModel;
drawGraphModel=function(){oldGraphDraw();graphicalResultsOverlay()};
const oldAnalyzeCad=$("cadAnalyze");
if(oldAnalyzeCad){oldAnalyzeCad.onclick=()=>analyzeGraphicalModel()}

function designGraphicalContinuousBeam(){
 try{
  analyzeGraphicalModel();
  const fc=num("bfc")||28,fy=num("bfy")||420,cover=num("bcover")||40,db=num("bbar")||16,ds=num("bstir")||8;
  const b=num("bb"),h=num("bh"),d=h-cover-ds-db/2,phi=.90,area=Math.PI*db*db/4;
  if(d<=0)throw new Error("Invalid beam effective depth.");
  const spans=(lastBeamResult&&lastBeamResult.spans)||[];
  const amin=Math.max(1.4/fy*b*d,.25*Math.sqrt(fc)/fy*b*d);
  const designs=spans.map((sp,i)=>{
   const pos=Math.max(0,Number(sp.maxPositiveMoment||sp.Mpos||0));
   const neg=Math.max(0,Number(sp.maxNegativeMoment||Math.abs(sp.Mneg||0)));
   function calc(M){
    const R=M*1e6/(phi*b*d*d),disc=1-2*R/(.85*fc);
    if(disc<0)return {n:2,As:2*area,over:true};
    const req=Math.max((.85*fc*b/fy)*(d-d*Math.sqrt(disc)),amin);
    const n=Math.max(2,Math.ceil(req/area));
    return {n,As:n*area,over:false};
   }
   return {span:i+1,bottom:calc(pos),top:calc(neg)};
  });
  window.lastContinuousDesign={fc,fy,b,h,d,db,ds,designs};
  const p=document.querySelector(".property-panel");
  if(p){let x=document.getElementById("continuousDesign");if(!x){x=document.createElement("div");x.id="continuousDesign";x.className="graph-summary";p.appendChild(x)}
   x.innerHTML="<b>PRELIMINARY RCC REINFORCEMENT</b>"+designs.map(v=>"<span>Span "+v.span+" • Bottom "+v.bottom.n+"Ø"+db+" • Top "+v.top.n+"Ø"+db+"</span>").join("");
  }
  setCadMode("design");drawReinforcementOverlay();
 }catch(e){alert(e.message)}
}
function drawReinforcementOverlay(){
 const c=$("beamModelCanvas");if(!c||!window.lastContinuousDesign)return;
 const ctx=c.getContext("2d"),g=graphCanvasGeometry(),nodes=cadNodes.slice().sort((a,b)=>a.x-b.x);
 ctx.save();ctx.lineWidth=3;ctx.font="bold 10px system-ui";
 for(let i=0;i<nodes.length-1;i++){
  const x1=g.x1+nodes[i].x/g.L*g.s,x2=g.x1+nodes[i+1].x/g.L*g.s;
  ctx.beginPath();ctx.moveTo(x1+12,g.y+10);ctx.lineTo(x2-12,g.y+10);ctx.stroke();
  ctx.fillText("BOTTOM +M",(x1+x2)/2-30,g.y+27);
 }
 for(let i=1;i<nodes.length-1;i++){
  const x=g.x1+nodes[i].x/g.L*g.s;
  ctx.beginPath();ctx.moveTo(x-24,g.y-9);ctx.lineTo(x+24,g.y-9);ctx.stroke();
  ctx.fillText("TOP −M",x-25,g.y-22);
 }
 ctx.restore();
}

function drawBeamDetailing(){
 const c=$("beamModelCanvas");if(!c||!window.lastContinuousDesign)return;
 const ctx=c.getContext("2d"),g=graphCanvasGeometry(),nodes=cadNodes.slice().sort((a,b)=>a.x-b.x),D=window.lastContinuousDesign,y=g.y;
 ctx.save();ctx.lineWidth=2;ctx.font="10px system-ui";
 for(let i=0;i<nodes.length-1;i++){const x1=g.x1+nodes[i].x/g.L*g.s,x2=g.x1+nodes[i+1].x/g.L*g.s,bd=D.designs[i];
  ctx.beginPath();ctx.moveTo(x1+10,y+5);ctx.lineTo(x2-10,y+5);ctx.stroke();ctx.fillText("B"+bd.bottom.n+"Ø"+D.db,(x1+x2)/2-22,y+34);
  ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(x1+8,y-16);ctx.lineTo(Math.min(x1+60,x2-20),y-16);ctx.stroke();ctx.beginPath();ctx.moveTo(Math.max(x2-60,x1+20),y-16);ctx.lineTo(x2-8,y-16);ctx.stroke();ctx.setLineDash([]);
  ctx.fillText("STIRRUP ZONE",x1+12,y-22);
 }
 for(let i=1;i<nodes.length-1;i++){const x=g.x1+nodes[i].x/g.L*g.s,bd=D.designs[i-1];ctx.beginPath();ctx.moveTo(x-45,y-5);ctx.lineTo(x+45,y-5);ctx.stroke();ctx.fillText("T"+bd.top.n+"Ø"+D.db,x-25,y-19);}
 ctx.restore();
}
function showDetailingSchedule(){
 if(!window.lastContinuousDesign)return;const D=window.lastContinuousDesign;
 let p=document.getElementById("detailingSchedule");
 if(!p){p=document.createElement("div");p.id="detailingSchedule";p.className="graph-summary";document.querySelector(".property-panel")?.appendChild(p)}
 p.innerHTML="<b>BAR SCHEDULE • PRELIMINARY</b>"+D.designs.map(x=>"<span>Span "+x.span+": B"+x.bottom.n+"Ø"+D.db+" bottom; top over internal support; closer stirrup zones near supports</span>").join("");
}
const _drawGraphModelBase=drawGraphModel;
drawGraphModel=function(){_drawGraphModelBase();if(window.lastContinuousDesign){drawBeamDetailing();showDetailingSchedule();}};

/* Beam action-point workflow — functional controls */
function initBeamActionPoints(){
 const root=document.getElementById("beamActionPoints"); if(!root)return;
 const loadsBody=document.getElementById("actionLoadsBody");
 const stiffBody=document.getElementById("actionStiffnessBody");
 const contBody=document.getElementById("actionContinuousBody");
 const designBody=document.getElementById("actionDesignBody");
 const makeBtn=(text,cls,fn)=>{const b=document.createElement("button");b.className=cls||"calculate";b.textContent=text;b.type="button";b.onclick=fn;return b};
 const note=(txt)=>{const p=document.createElement("p");p.className="action-note";p.textContent=txt;return p};

 // Move the real analysis/design cards into their action panels.
 const beam=document.getElementById("beam");
 const grid=beam?.querySelector(".grid");
 const input=grid?.querySelector(".card.inputs"), result=grid?.querySelector(".card.results");
 if(input)stiffBody.appendChild(input);
 if(result)stiffBody.appendChild(result);
 const continuous=beam?.querySelector(".continuous-beam");
 if(continuous)contBody.appendChild(continuous);
 const cbNext=continuous?.nextElementSibling;
 if(cbNext && cbNext.classList.contains("matrix-card"))contBody.appendChild(cbNext);
 const diagram=document.getElementById("cbDiagram");
 if(diagram)contBody.appendChild(diagram);
 const design=beam?.querySelector(".beam-design");
 if(design)designBody.appendChild(design);
 const matrix=[...beam.querySelectorAll(".matrix-card")];
 const singleMatrix=matrix.find(x=>!x.classList.contains("continuous-beam")&&!x.classList.contains("beam-design")&&x!==cbNext);
 if(singleMatrix)stiffBody.appendChild(singleMatrix);

 // 01 — actual load calculation, with results shown inline.
 const loadBox=document.createElement("div");loadBox.className="action-linked-tool";
 const loadText=document.createElement("div");loadText.innerHTML="<b>Calculate building loads</b><small>Runs the BNBC load calculator using the current Load Calculator inputs.</small>";
 const loadBtn=makeBtn("CALCULATE LOADS →","action-open-load",()=>{calculate();document.getElementById("loads")?.scrollIntoView({behavior:"smooth",block:"start"});});
 loadBox.append(loadText,loadBtn);loadsBody.appendChild(loadBox);
 loadsBody.appendChild(note("The button runs the real load-calculation function. The detailed D, L and factored results remain in Load Calculator."));

 // 02 — actual one-span stiffness run.
 const stiffBox=document.createElement("div");stiffBox.className="action-linked-tool";
 const stiffText=document.createElement("div");stiffText.innerHTML="<b>Run direct stiffness analysis</b><small>Solves [K]{Δ}={P}, reactions, moments, shear and displacement.</small>";
 const stiffBtn=makeBtn("RUN STIFFNESS →","action-open-load",()=>{beamStiffness();updateCad();document.getElementById("actionStiffness")?.setAttribute("open","");});
 stiffBox.append(stiffText,stiffBtn);stiffBody.insertBefore(stiffBox,stiffBody.firstChild);

 // 03 — actual multi-span continuous stiffness run.
 const contBox=document.createElement("div");contBox.className="action-linked-tool";
 const contText=document.createElement("div");contText.innerHTML="<b>Run continuous-beam stiffness</b><small>Solves the global multi-span system and updates span/joint actions and SFD/BMD.</small>";
 const contBtn=makeBtn("RUN CONTINUOUS ANALYSIS →","action-open-load",()=>{continuousBeamStiffness();document.getElementById("actionContinuous")?.setAttribute("open","");});
 contBox.append(contText,contBtn);contBody.insertBefore(contBox,contBody.firstChild);

 // 04 — actual reinforcement design.
 const designBox=document.createElement("div");designBox.className="action-linked-tool";
 const designText=document.createElement("div");designText.innerHTML="<b>Design reinforcement</b><small>Runs the existing preliminary RCC flexure/shear design from the solved beam actions.</small>";
 const designBtn=makeBtn("DESIGN REINFORCEMENT →","action-open-load",()=>{beamDesign();document.getElementById("actionDesign")?.setAttribute("open","");});
 designBox.append(designText,designBtn);designBody.insertBefore(designBox,designBody.firstChild);

 root.querySelectorAll("details").forEach(d=>d.addEventListener("toggle",()=>{if(d.open)d.scrollIntoView({behavior:"smooth",block:"nearest"})}));
}
initBeamActionPoints();

/* Unified app navigation: sidebar, dashboard cards and mobile tab bar */
(function initAppNavigation(){
 function go(section){
  document.querySelectorAll(".section").forEach(x=>x.classList.toggle("active-section",x.id===section));
  document.querySelectorAll("[data-section]").forEach(x=>x.classList.toggle("active",x.dataset.section===section));
  const title=document.querySelector(".topbar h1");
  const sub=document.querySelector(".topbar .subtitle");
  const meta={home:["Structural Design Toolkit","Loads, analysis and RCC design in one workspace."],loads:["Building Load Calculator","Slab, beam, wall, finishing and occupancy-based load calculations."],beam:["RCC Beam Workbench","Stiffness analysis, continuous beam actions and preliminary reinforcement."],column:["RCC Column Design","Axial load, interaction and preliminary reinforcement design."],slab:["RCC Slab Design","One-way and two-way slab analysis and preliminary reinforcement."]}[section];
  if(meta&&title){title.textContent=meta[0];if(sub)sub.textContent=meta[1]}
  window.scrollTo({top:0,behavior:"smooth"});
 }
 document.addEventListener("click",e=>{const el=e.target.closest("[data-section]");if(el){e.preventDefault();go(el.dataset.section)}});
 go(document.querySelector(".section.active-section")?.id||"home");
})();

(function initColumnPreview(){const b=$("columnPreview");if(!b)return;b.onclick=()=>{const w=num("cbw"),h=num("cbh"),fc=num("cfc"),fy=num("cfy"),Pu=num("cpu"),Mu=num("cmu"),Ag=w*h,stress=Pu*1000/Ag,bar=num("cbar"),Ast=8*Math.PI*bar*bar/4,rho=Ast/Ag*100,ok=stress<=0.35*fc;$("cAg").textContent=Ag.toLocaleString();$("cPuOut").textContent=Pu.toFixed(1);$("cMuOut").textContent=Mu.toFixed(1);$("cStress").textContent=stress.toFixed(2);$("cRho").textContent=rho.toFixed(2);$("cStatus").textContent=ok?"TRIAL OK — continue design":"REVIEW SECTION";};})();

(function initColumnInteraction(){
 const run=$("columnInteraction");if(!run)return;
 run.onclick=function(){
  const w=num("cbw"),h=num("cbh"),fc=num("cfc"),fy=num("cfy"),Pu=num("cpu"),Mu=num("cmu"),Ag=w*h,db=num("cbar"),Ast=8*Math.PI*db*db/4,phi=.65,Po=.85*fc*(Ag-Ast)+fy*Ast,points=[];
  for(let i=0;i<=20;i++){const r=i/20;points.push([Mu*1.8*(1-Math.pow(r,.7)),phi*Po*(1-r*.72)]);}
  const c=$("columnInteractionCanvas"),ctx=c.getContext("2d"),W=c.width,H=c.height,maxM=Math.max.apply(null,points.map(function(p){return p[0]}),Mu*1.2),maxP=Math.max.apply(null,points.map(function(p){return p[1]}),Pu*1.2),X=function(m){return 45+m/maxM*(W-75)},Y=function(p){return H-35-p/maxP*(H-65)};
  ctx.clearRect(0,0,W,H);ctx.beginPath();points.forEach(function(p,i){if(i)ctx.lineTo(X(p[0]),Y(p[1]));else ctx.moveTo(X(p[0]),Y(p[1]))});ctx.strokeStyle="#007aff";ctx.lineWidth=3;ctx.stroke();ctx.beginPath();ctx.moveTo(X(0),Y(0));ctx.lineTo(X(Mu),Y(Pu));ctx.strokeStyle="#ff3b30";ctx.setLineDash([6,5]);ctx.stroke();ctx.setLineDash([]);ctx.fillStyle="#111827";ctx.font="12px system-ui";ctx.fillText("M (kN·m)",W-90,H-10);ctx.fillText("P (kN)",8,22);
  $("interactionStatus").textContent=Pu<=phi*Po*.28?"Trial point inside curve":"Review section";$("interactionSummary").textContent="Preliminary screening curve. Final ACI/BNBC interaction analysis requires strain compatibility, phi transition and all applicable load combinations.";
 };
 const bars=$("columnBars");if(bars)bars.onclick=function(){const Ag=num("cbw")*num("cbh"),opts=[12,16,20,25,28],min=.01*Ag,max=.08*Ag,out=[];opts.forEach(function(d){const A=12*Math.PI*d*d/4;if(A>=min&&A<=max)out.push("12Ø"+d+" • "+A.toFixed(0)+" mm²")});$("columnBarOptions").textContent=out.length?out.join("   |   "):"No trial bar option in the preliminary 1–8% range.";};
})();
(function initTwoWay(){const b=$("twoWayAnalyze");if(!b)return;b.onclick=function(){const Lx=num("twLx"),Ly=num("twLy"),D=num("twd"),L=num("twl"),wu=Math.max(1.4*D,1.2*D+1.6*L),ratio=Lx/Ly,Mo=wu*Ly*Ly*Lx/8;$("twratio").textContent=ratio.toFixed(2);$("twwu").textContent=wu.toFixed(2);$("twMo").textContent=Mo.toFixed(2);$("twshort").textContent=ratio<=2?"Strong two-way action":"Long-span direction becomes dominant";};})();
