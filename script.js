const $=id=>document.getElementById(id);function num(id){return Number($(id).value)||0}
function calculate(){
 const t=num("slab_t")/1000,gc=num("conc_w"),finish=num("finish"),live=num("live");
 const bb=num("beam_b")/1000,bh=num("beam_h")/1000,bl=num("beam_L");
 const wt=num("wall_t")/1000,wh=num("wall_h"),ww=num("wall_w"),pl=num("plaster");
 const slab=t*gc,dead=slab+finish,service=dead+live;
 const beam=bb*bh*gc;
 const wall=wt*wh*ww+2*wh*pl;
 const factored=Math.max(1.4*dead,1.2*dead+1.6*live);
 $("slab_out").textContent=slab.toFixed(2);$("finish_out").textContent=finish.toFixed(2);
 $("dead_out").textContent=dead.toFixed(2);$("live_out").textContent=live.toFixed(2);
 $("service_out").textContent=service.toFixed(2);$("beam_out").textContent=beam.toFixed(2);
 $("wall_out").textContent=wall.toFixed(2);$("factored_out").textContent=factored.toFixed(2);
}
$("calculate").addEventListener("click",calculate);document.querySelectorAll("input").forEach(i=>i.addEventListener("input",calculate));
document.querySelectorAll(".nav-item").forEach(btn=>btn.addEventListener("click",()=>{document.querySelectorAll(".nav-item").forEach(x=>x.classList.remove("active"));btn.classList.add("active");document.querySelectorAll(".section,.placeholder").forEach(x=>x.style.display="none");const t=$(btn.dataset.section);if(t)t.style.display="block";}));
calculate();