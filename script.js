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