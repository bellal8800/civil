const $=id=>document.getElementById(id);
function num(id){return Number($(id).value)||0}
function calculate(){
  const L=num("L"), b=num("b"), h=num("h"), cover=num("cover"), db=num("db"), ds=num("ds"), Dadd=num("Dadd"), LL=num("LL");
  if(L<=0||b<=0||h<=0){alert("Please enter valid beam dimensions.");return}
  const sw=b*h*24;
  const D=sw+Dadd;
  const wu=Math.max(1.4*D,1.2*D+1.6*LL);
  const Mu=wu*L*L/8;
  const Vu=wu*L/2;
  const d=h*1000-cover-ds-db/2;
  $("sw").textContent=sw.toFixed(2);
  $("D").textContent=D.toFixed(2);
  $("wu").textContent=wu.toFixed(2);
  $("Mu").textContent=Mu.toFixed(2);
  $("Vu").textContent=Vu.toFixed(2);
  $("d").textContent=Math.max(d,0).toFixed(0);
}
$("calculate").addEventListener("click",calculate);
document.querySelectorAll(".field input").forEach(i=>i.addEventListener("input",calculate));
document.querySelectorAll(".nav-item").forEach(btn=>{
  btn.addEventListener("click",()=>{
    document.querySelectorAll(".nav-item").forEach(x=>x.classList.remove("active"));
    btn.classList.add("active");
    document.querySelectorAll(".section,.placeholder").forEach(x=>x.style.display="none");
    const target=$(btn.dataset.section);
    if(target) target.style.display="block";
  });
});
calculate();