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
$("calculate").addEventListener("click",calculate);
$("slabCalculate").addEventListener("click",slabDesign);
document.querySelectorAll("#slab input,#slab select").forEach(i=>i.addEventListener("input",slabDesign));
document.querySelectorAll("input,select").forEach(i=>i.addEventListener("input",calculate));
document.querySelectorAll(".nav-item").forEach(btn=>btn.addEventListener("click",()=>{
 document.querySelectorAll(".nav-item").forEach(x=>x.classList.remove("active"));btn.classList.add("active");
 document.querySelectorAll(".section,.placeholder").forEach(x=>x.style.display="none");
 const t=$(btn.dataset.section);if(t)t.style.display="block";
}));
calculate();