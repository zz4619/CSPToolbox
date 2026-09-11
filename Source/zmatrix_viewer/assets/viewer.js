// Bundled into each export. All state and controls belong to this viewer root.
const root = document.getElementById(rootId);
if (!root || root.dataset.ready === "true") return;
const get = role => root.querySelector(`[data-role="${role}"]`);
const fail = error => {
  get("loading").hidden = true;
  const message = get("error");
  message.hidden = false;
  message.textContent = `Viewer could not render: ${error.message || error}`;
  root.dataset.ready = "error";
};
try {
  const data = JSON.parse(get("data").textContent);
  if (data.schema_version !== 1 || !data.scenes.length) throw new Error("Unsupported or empty viewer data.");
  const canvas = get("canvas"), ctx = canvas.getContext("2d");
  if (!ctx) throw new Error("This browser does not provide a 2D canvas.");
  const state = {scene: 0, selected: null, hovered: null, hoveredAtom: null, rotX: 0, rotY: 0,
    zoom: 1, dragging: false, lastX: 0, lastY: 0, center: [0,0,0], radius: 1, renderCount: 0};
  const listeners = [];
  const on = (node, name, fn, options) => {
    const guarded = event => {try {fn(event);} catch (error) {fail(error);}};
    node.addEventListener(name, guarded, options);
    listeners.push(() => node.removeEventListener(name, guarded, options));
  };
  const scene = () => data.scenes[state.scene];
  const molecule = () => get("cell").value === "full" && scene().full_cell ? scene().full_cell : scene().molecule;
  const active = () => scene().coordinates.find(c => c.name === (state.hovered || state.selected));
  const shownCoordinates = () => scene().coordinates.filter(c => get("filter").value === "all" ||
    (get("filter").value === "independent" ? c.independent : c.kind === "dih"));
  const linkKey = (a,b) => a<b ? `${a}-${b}` : `${b}-${a}`;
  const delta = c => c.kind === "dih" ? (((c.value-c.reference_value+180)%360)+360)%360-180 : c.value-c.reference_value;
  const number = (value, unit) => `${value.toFixed(unit === "Å" ? 4 : 2)}${unit}`;
  const colorProbe = document.createElement("span");
  colorProbe.hidden = true;
  root.appendChild(colorProbe);
  let palette;
  function colors() {
    if (palette) return palette;
    const resolve = token => {
      colorProbe.style.color = `var(--csp-${token})`;
      return getComputedStyle(colorProbe).color;
    };
    palette = {bg:resolve("bg"),fg:resolve("fg"),muted:resolve("muted"),border:resolve("border"),
      chains:[1,2,3,4,5].map(i=>resolve(`torsion-${i}`))};
    return palette;
  }
  function chainColor(coordinate) {
    const independent = scene().coordinates.filter(c=>c.independent);
    const index = independent.findIndex(c=>c.name === coordinate?.name);
    return colors().chains[Math.max(0,index)%5];
  }
  function addText(parent, tag, text, className) {
    const element = document.createElement(tag);
    element.textContent = text;
    if (className) element.className = className;
    parent.appendChild(element);
    return element;
  }
  function buildList() {
    get("coordinates").replaceChildren();
    for (const c of shownCoordinates()) {
      const button = document.createElement("button");
      button.type = "button"; button.className = "csp-coordinate"; button.dataset.coordinate = c.name;
      addText(button,"span",c.name,"csp-name");
      addText(button,"span",number(c.value,c.unit),"csp-value");
      addText(button,"small",c.atom_labels.join("–") + (c.independent ? " · independent" : ""));
      if (c.reference_value !== null) addText(button,"small",`Reference ${number(c.reference_value,c.unit)} · Δ ${number(delta(c),c.unit)}`);
      get("coordinates").appendChild(button);
    }
    if (!shownCoordinates().length) addText(get("coordinates"),"p",scene().coordinates.length ?
      "No coordinates in this selection." : "No Z-matrix supplied. The molecular geometry is displayed directly.");
    updateSelection();
  }
  function updateSelection() {
    const coordinate = active();
    for (const button of get("coordinates").querySelectorAll("button")) {
      button.setAttribute("aria-pressed", String(button.dataset.coordinate === state.selected));
    }
    get("detail").replaceChildren();
    if (!coordinate) return;
    addText(get("detail"),"strong",`${coordinate.name}: ${coordinate.atom_labels.join("–")}`);
    addText(get("detail"),"div",number(coordinate.value,coordinate.unit));
    if (coordinate.reference_value !== null) addText(get("detail"),"div",`Reference ${number(coordinate.reference_value,coordinate.unit)} · Δ ${number(delta(coordinate),coordinate.unit)}`);
    if (coordinate.lower !== null && coordinate.upper !== null) addText(get("detail"),"div",`Defined interval: ${number(coordinate.lower,coordinate.unit)} to ${number(coordinate.upper,coordinate.unit)}`);
  }
  function setupScene() {
    const s = scene();
    get("title").textContent = s.title;
    get("subtitle").textContent = s.mapping_method || (s.full_cell ? "Crystal structure" : "Molecular geometry");
    get("reference-control").hidden = !s.reference;
    get("cell-control").hidden = !s.full_cell;
    get("cell").value = "asu"; get("reference").checked = true;
    const hasIndependent = s.coordinates.some(c=>c.independent);
    get("filter").querySelector('[value="independent"]').disabled = !hasIndependent;
    get("filter").value = hasIndependent ? "independent" : "dihedrals";
    get("filter").disabled = !s.coordinates.length;
    state.selected = shownCoordinates()[0]?.name || null; state.hovered = null;
    get("mapping-control").hidden = !s.mapping.length;
    get("mapping").replaceChildren();
    if (s.mapping.length) {
      const table = document.createElement("table");
      const header = document.createElement("tr");
      addText(header,"th","Reference / Z-matrix"); addText(header,"th","Input atom"); table.appendChild(header);
      for (const [target,source] of s.mapping) {
        const row = document.createElement("tr"); addText(row,"td",target); addText(row,"td",source); table.appendChild(row);
      }
      get("mapping").appendChild(table);
    }
    root.querySelector('[data-action="download-zmatrix"]').hidden = !s.mapped_zmatrix;
    get("notices").textContent = [...s.notices,...s.molecule.warnings].join(" ");
    get("legend").textContent = s.reference ? "Colour: comparison · Grey: reference · Dashed: internal-coordinate reference link" :
      (s.coordinates.length ? "Selected coordinate: highlighted atom chain · Dashed: reference link, not a covalent bond" : "Bonds inferred by distance unless supplied in the input");
    buildList(); resetView();
  }
  function fit() {
    const points = [...molecule().atoms];
    if (scene().reference) points.push(...scene().reference.atoms);
    if (scene().cell.length) points.push(...cellCorners());
    state.center = ["x","y","z"].map(key=>points.reduce((sum,p)=>sum+p[key],0)/points.length);
    state.radius = Math.max(1,...points.map(p=>Math.hypot(p.x-state.center[0],p.y-state.center[1],p.z-state.center[2])));
  }
  function resetView() {
    state.zoom = 1; fit();
    const variance = ["x","y","z"].map((key,i)=>molecule().atoms.reduce((sum,p)=>sum+(p[key]-state.center[i])**2,0));
    const small = variance.indexOf(Math.min(...variance));
    state.rotX = small === 1 ? 1.25 : 0.15;
    state.rotY = small === 0 ? 1.25 : 0.15;
    state.hoveredAtom = null; resize();
  }
  function resize() {
    const rect = canvas.getBoundingClientRect(), ratio = window.devicePixelRatio || 1;
    canvas.width = Math.max(1,Math.round(rect.width*ratio)); canvas.height = Math.max(1,Math.round(rect.height*ratio));
    ctx.setTransform(ratio,0,0,ratio,0,0); draw();
  }
  const scale = () => Math.min(canvas.clientWidth,canvas.clientHeight)/(state.radius*2.65)*state.zoom;
  function project(p) {
    const x=p.x-state.center[0],y=p.y-state.center[1],z=p.z-state.center[2];
    const cy=Math.cos(state.rotY),sy=Math.sin(state.rotY),cx=Math.cos(state.rotX),sx=Math.sin(state.rotX);
    const xx=x*cy+z*sy, zz=-x*sy+z*cy;
    return {x:canvas.clientWidth/2+xx*scale(),y:canvas.clientHeight/2-(y*cx-zz*sx)*scale(),z:y*sx+zz*cx};
  }
  function cellCorners() {
    return Array.from({length:8},(_,mask)=>{
      const p=[0,0,0]; for(let axis=0;axis<3;axis++) if(mask&(1<<axis)) for(let k=0;k<3;k++) p[k]+=scene().cell[axis][k];
      return {x:p[0],y:p[1],z:p[2]};
    });
  }
  function line(a,b,color,width,dashed=false,alpha=1) {
    ctx.save();ctx.globalAlpha=alpha;ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);
    ctx.lineWidth=width;ctx.lineCap="round";ctx.strokeStyle=color;if(dashed)ctx.setLineDash([5,4]);ctx.stroke();ctx.restore();
  }
  function drawMolecule(mol, ghost) {
    const c=active(), activeIds=new Set(!ghost && c ? c.atom_indices : []);
    const links=new Set(); if(c && !ghost) for(let i=1;i<c.atom_indices.length;i++) links.add(linkKey(c.atom_indices[i-1],c.atom_indices[i]));
    const atoms=mol.atoms.filter(a=>get("hydrogens").checked || a.element!=="H" || activeIds.has(a.index));
    const projected=new Map(atoms.map(a=>[a.index,project(a)]));
    const chain=chainColor(c), color=colors();
    const bonds=mol.bonds.filter(b=>projected.has(b.left)&&projected.has(b.right));
    bonds.sort((a,b)=>(projected.get(a.left).z+projected.get(a.right).z)-(projected.get(b.left).z+projected.get(b.right).z));
    for(const b of bonds) line(projected.get(b.left),projected.get(b.right),ghost?color.muted:links.has(linkKey(b.left,b.right))?chain:color.muted,
      ghost?2.8:links.has(linkKey(b.left,b.right))?5:3.2,false,ghost?0.28:1);
    if(c&&!ghost) {
      const bonded=new Set(bonds.map(b=>linkKey(b.left,b.right)));
      for(let i=1;i<c.atom_indices.length;i++) {
        const a=c.atom_indices[i-1],b=c.atom_indices[i];
        if(!bonded.has(linkKey(a,b))&&projected.has(a)&&projected.has(b)) line(projected.get(a),projected.get(b),chain,2.5,true);
      }
    }
    const visible=atoms.map(atom=>({atom,p:projected.get(atom.index)})).sort((a,b)=>a.p.z-b.p.z);
    for(const {atom,p} of visible) {
      const r=Math.max(3.5,Math.min(15,atom.display_radius*scale()*.8));
      ctx.save();ctx.globalAlpha=ghost?0.25:1;
      if(activeIds.has(atom.index)||(!ghost&&state.hoveredAtom===atom.index)) {
        ctx.beginPath();ctx.arc(p.x,p.y,r+4,0,Math.PI*2);ctx.strokeStyle=chain;ctx.lineWidth=2;ctx.stroke();
      }
      ctx.beginPath();ctx.arc(p.x,p.y,r,0,Math.PI*2);
      ctx.fillStyle=ghost?color.muted:atom.element==='C'?color.fg:atom.element==='H'?color.bg:atom.color;
      ctx.fill();ctx.strokeStyle=color.border;ctx.lineWidth=1;ctx.stroke();ctx.restore();
    }
    if(!ghost&&get("labels").checked) drawLabels(visible,activeIds);
  }
  function drawLabels(visible,activeIds) {
    const boxes=[],col=colors();ctx.font="12px system-ui,sans-serif";ctx.textBaseline="middle";
    const overlap=(a,b)=>Math.max(0,Math.min(a.x+a.w,b.x+b.w)-Math.max(a.x,b.x))*Math.max(0,Math.min(a.y+a.h,b.y+b.h)-Math.max(a.y,b.y));
    for(const {atom,p} of [...visible].sort((a,b)=>Number(activeIds.has(b.atom.index))-Number(activeIds.has(a.atom.index))||a.p.y-b.p.y)) {
      if(visible.length>160&&!activeIds.has(atom.index)&&state.hoveredAtom!==atom.index) continue;
      const w=ctx.measureText(atom.label).width+4;let best;
      for(const pad of [14,23,32]) for(const [dx,dy] of [[1,-1],[-1,-1],[1,1],[-1,1],[1,0],[-1,0],[0,-1],[0,1]]) {
        const b={x:p.x+dx*pad-(dx<0?w:dx===0?w/2:0),y:p.y+dy*pad-7,w,h:14};let score=pad*.02;
        if(b.x<2||b.y<2||b.x+w>canvas.clientWidth-2||b.y+14>canvas.clientHeight-2)score+=10000;
        for(const old of boxes)score+=overlap(b,old)*10;
        for(const other of visible) score+=overlap(b,{x:other.p.x-8,y:other.p.y-8,w:16,h:16});
        if(!best||score<best.score)best={...b,score};
      }
      boxes.push(best);ctx.lineWidth=3;ctx.strokeStyle=col.bg;ctx.strokeText(atom.label,best.x+2,best.y+7);
      ctx.fillStyle=col.fg;ctx.fillText(atom.label,best.x+2,best.y+7);
    }
  }
  function draw() {
    if(!canvas.clientWidth||!canvas.clientHeight)return;
    ctx.clearRect(0,0,canvas.clientWidth,canvas.clientHeight);
    if(scene().cell.length) {
      const corners=cellCorners().map(project);
      for(let i=0;i<8;i++)for(let a=0;a<3;a++){const j=i^(1<<a);if(i<j)line(corners[i],corners[j],colors().border,1.2);}
    }
    if(scene().reference&&get("reference").checked)drawMolecule(scene().reference,true);
    drawMolecule(molecule(),false);state.renderCount++;
    canvas.setAttribute("aria-description",`View rotated ${(state.rotX*180/Math.PI).toFixed(1)} and ${(state.rotY*180/Math.PI).toFixed(1)} degrees; zoom ${(state.zoom*100).toFixed(0)} percent.`);
    get("status").textContent=`${molecule().atoms.length} atoms · ${molecule().bonds.length} bonds`+
      (scene().rmsd===null?"":` · Aligned all-atom RMSD ${scene().rmsd.toFixed(4)} Å`)+
      (molecule().atoms.length>160?" · Hover to label individual atoms":"");
  }
  function hitTestAtom(event) {
    const rect=canvas.getBoundingClientRect();let best=null,distance=Infinity;
    for(const atom of molecule().atoms) {
      if(atom.element==='H'&&!get("hydrogens").checked)continue;
      const p=project(atom),d=Math.hypot(p.x-(event.clientX-rect.left),p.y-(event.clientY-rect.top));
      if(d<14&&d<distance){best=atom.index;distance=d;}
    }
    return best;
  }
  on(root,"click",event=>{
    const coordinate=event.target.closest("[data-coordinate]");
    if(coordinate&&root.contains(coordinate)) {
      state.selected=state.selected===coordinate.dataset.coordinate?null:coordinate.dataset.coordinate;
      state.hovered=null;updateSelection();draw();return;
    }
    const action=event.target.closest("[data-action]")?.dataset.action;
    if(action==="reset")resetView();
    if(action==="clear"){state.selected=null;state.hovered=null;state.hoveredAtom=null;updateSelection();draw();}
    if(action==="download-zmatrix"&&scene().mapped_zmatrix) {
      const url=URL.createObjectURL(new Blob([scene().mapped_zmatrix],{type:"text/plain"}));
      const a=document.createElement("a");a.href=url;a.download="mapped-molecule.zmat";a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
    }
  });
  on(get("filter"),"change",()=>{state.hovered=null;if(!shownCoordinates().some(c=>c.name===state.selected))state.selected=null;buildList();draw();});
  on(get("scene"),"change",()=>{state.scene=Number(get("scene").value);setupScene();});
  on(get("cell"),"change",resetView);
  for(const role of ["reference","hydrogens","labels"])on(get(role),"change",draw);
  on(root,"keydown",event=>{if(event.key==="Escape"){state.selected=null;state.hovered=null;updateSelection();draw();}});
  on(canvas,"keydown",event=>{
    if(!["ArrowLeft","ArrowRight","ArrowUp","ArrowDown","+","=","-","Home"].includes(event.key))return;
    event.preventDefault();
    if(event.key==="Home"){resetView();return;}
    if(event.key==="ArrowLeft")state.rotY-=.1;
    if(event.key==="ArrowRight")state.rotY+=.1;
    if(event.key==="ArrowUp")state.rotX-=.1;
    if(event.key==="ArrowDown")state.rotX+=.1;
    if(["+","=","-"].includes(event.key))state.zoom=Math.max(.35,Math.min(5,state.zoom*(event.key==="-"?1/1.1:1.1)));
    draw();
  });
  on(canvas,"pointerdown",event=>{state.dragging=true;state.lastX=event.clientX;state.lastY=event.clientY;canvas.setPointerCapture(event.pointerId);canvas.classList.add("is-dragging");});
  on(canvas,"pointermove",event=>{
    if(state.dragging){state.rotY+=(event.clientX-state.lastX)*.01;state.rotX+=(event.clientY-state.lastY)*.01;state.lastX=event.clientX;state.lastY=event.clientY;}
    else{state.hoveredAtom=hitTestAtom(event);canvas.classList.toggle("is-atom-hover",state.hoveredAtom!==null);}
    draw();
  });
  const stop=event=>{state.dragging=false;canvas.classList.remove("is-dragging");if(canvas.hasPointerCapture(event.pointerId))canvas.releasePointerCapture(event.pointerId);};
  on(canvas,"pointerup",stop);on(canvas,"pointercancel",stop);
  on(canvas,"pointerleave",()=>{if(!state.dragging){state.hoveredAtom=null;draw();}});
  on(canvas,"wheel",event=>{event.preventDefault();state.zoom=Math.max(.35,Math.min(5,state.zoom*Math.exp(-event.deltaY*.001)));draw();},{passive:false});
  const observer=new ResizeObserver(()=>{try{resize();}catch(error){fail(error);}});observer.observe(canvas);
  const themeObserver=new MutationObserver(()=>{palette=null;draw();});themeObserver.observe(document.documentElement,{attributes:true,attributeFilter:["style","class","data-theme"]});
  const media=window.matchMedia("(prefers-color-scheme: dark)");on(media,"change",()=>{palette=null;draw();});
  for(let i=0;i<data.scenes.length;i++){const option=document.createElement("option");option.value=String(i);option.textContent=data.scenes[i].title;get("scene").appendChild(option);}
  get("scene-control").hidden=data.scenes.length<2;
  root.cspViewer={getState:()=>({scene:state.scene,selected:state.selected,rotation:[state.rotX,state.rotY],zoom:state.zoom,
    cell:get("cell").value,referenceVisible:get("reference").checked,renderCount:state.renderCount}),
    destroy:()=>{listeners.forEach(remove=>remove());observer.disconnect();themeObserver.disconnect();colorProbe.remove();delete root.cspViewer;delete root.dataset.ready;}};
  setupScene();get("loading").hidden=true;root.dataset.ready="true";
} catch(error) {fail(error);}
