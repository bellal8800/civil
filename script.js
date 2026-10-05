const KEY="civilWorkspaceV1";
let state=JSON.parse(localStorage.getItem(KEY)||"null")||{folders:[
{id:"beam",name:"RCC Beam",open:true,tools:[{id:"beam-analysis",name:"Beam Analysis"}]},
{id:"column",name:"RCC Column",open:false,tools:[]},
{id:"slab",name:"RCC Slab",open:false,tools:[]},
{id:"foundation",name:"Foundation",open:false,tools:[]},
{id:"other",name:"Other Tools",open:false,tools:[]}]};
let selected=null,mode=null,parentId=null;
const tree=document.querySelector("#tree"),modal=document.querySelector("#modal"),input=document.querySelector("#modalInput");
function save(){localStorage.setItem(KEY,JSON.stringify(state))}
function esc(s){return s.replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]))}
function render(){
 const q=document.querySelector("#search").value.toLowerCase(); tree.innerHTML="";
 state.folders.forEach(f=>{
  if(q&&!((f.name+" "+f.tools.map(t=>t.name).join(" ")).toLowerCase().includes(q)))return;
  const wrap=document.createElement("div"),row=document.createElement("div"); row.className="row"+(selected===f.id?" selected":"");
  row.innerHTML='<span class="caret">'+(f.open?"▾":"▸")+'</span><span class="icon">📁</span><span class="name">'+esc(f.name)+'</span><span class="count">'+f.tools.length+'</span><button class="more">⋯</button>';
  row.onclick=e=>{if(e.target.classList.contains("more"))return;f.open=!f.open;selected=f.id;showFolder(f)};
  row.querySelector(".more").onclick=e=>{e.stopPropagation();folderMenu(f)}; wrap.appendChild(row);
  if(f.open){const ch=document.createElement("div");ch.className="children";
   f.tools.filter(t=>!q||t.name.toLowerCase().includes(q)).forEach(t=>{
    const tr=document.createElement("div");tr.className="row"+(selected===t.id?" selected":"");
    tr.innerHTML='<span class="caret"></span><span class="icon">🧮</span><span class="name">'+esc(t.name)+'</span>';
    tr.onclick=()=>showTool(f,t);ch.appendChild(tr)});wrap.appendChild(ch)}
  tree.appendChild(wrap)
 })
}
function showFolder(f){document.querySelector("#welcome").classList.add("hidden");document.querySelector("#toolPanel").classList.add("hidden");document.querySelector("#pageTitle").textContent=f.name;document.querySelector("#pageSub").textContent=f.tools.length+" tool"+(f.tools.length===1?"":"s")+" in this folder."}
function showTool(f,t){selected=t.id;document.querySelector("#welcome").classList.add("hidden");document.querySelector("#toolPanel").classList.remove("hidden");document.querySelector("#toolName").textContent=t.name;document.querySelector("#toolInfo").textContent="Folder: "+f.name;document.querySelector("#pageTitle").textContent=t.name;document.querySelector("#pageSub").textContent="Tool workspace • ACI 318-25";render()}
function openModal(m,pid=null){mode=m;parentId=pid;document.querySelector("#modalTitle").textContent=m==="folder"?"Create New Folder":"Create New Tool";document.querySelector("#modalHint").textContent=m==="folder"?"Enter the folder name.":"Enter the tool name.";document.querySelector("#saveModal").textContent=m==="folder"?"Create":"Add Tool";input.value="";modal.classList.remove("hidden");setTimeout(()=>input.focus(),20)}
function closeModal(){modal.classList.add("hidden")}
function submit(){const n=input.value.trim();if(!n)return;if(mode==="folder")state.folders.push({id:"f-"+Date.now(),name:n,open:true,tools:[]});else{const f=state.folders.find(x=>x.id===parentId);if(f)f.tools.push({id:"t-"+Date.now(),name:n})}save();closeModal();render()}
function folderMenu(f){const a=prompt("Folder action:
1 = Rename
2 = Delete
3 = Add Tool","1");if(a==="1"){const n=prompt("New folder name:",f.name);if(n&&n.trim()){f.name=n.trim();save();render()}}else if(a==="2"){if(confirm('Delete "'+f.name+'" and its tools?')){state.folders=state.folders.filter(x=>x.id!==f.id);save();render()}}else if(a==="3")openModal("tool",f.id)}
document.querySelector("#newFolderBtn").onclick=()=>openModal("folder");
document.querySelector("#welcomeFolder").onclick=()=>openModal("folder");
document.querySelector("#newToolBtn").onclick=()=>{const f=state.folders.find(x=>x.id===selected)||state.folders[0];openModal("tool",f.id)};
document.querySelector("#cancelModal").onclick=closeModal;document.querySelector("#saveModal").onclick=submit;
input.onkeydown=e=>{if(e.key==="Enter")submit();if(e.key==="Escape")closeModal()};
document.querySelector("#closeTool").onclick=()=>document.querySelector("#toolPanel").classList.add("hidden");
document.querySelector("#search").oninput=render;
document.querySelector("#resetBtn").onclick=()=>{if(confirm("Reset the workspace to the default folders?")){localStorage.removeItem(KEY);location.reload()}};
render();