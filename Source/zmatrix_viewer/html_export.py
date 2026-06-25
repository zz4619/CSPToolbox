"""Standalone HTML export for the Z-matrix viewer."""

from __future__ import annotations

from html import escape
import json
from pathlib import Path

from .model import ViewerMolecule


def write_viewer_html(molecule: ViewerMolecule, destination: str | Path) -> Path:
    """Write a standalone browser viewer for a reconstructed Z-matrix."""

    destination_path = Path(destination)
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    destination_path.write_text(render_viewer_html(molecule), encoding="utf-8")
    return destination_path


def render_viewer_html(molecule: ViewerMolecule) -> str:
    """Render a standalone HTML document containing viewer data and code."""

    data_json = json.dumps(_molecule_payload(molecule), ensure_ascii=True, allow_nan=False)
    data_json = data_json.replace("<", "\\u003c")
    title = escape(molecule.title or "Z-matrix viewer")
    return HTML_TEMPLATE.replace("__TITLE__", title).replace("__MOLECULE_DATA__", data_json)


def _molecule_payload(molecule: ViewerMolecule) -> dict[str, object]:
    return {
        "title": molecule.title,
        "source_name": molecule.source_name,
        "atoms": [
            {
                "index": atom.index,
                "label": atom.label,
                "element": atom.element,
                "x": atom.coordinates[0],
                "y": atom.coordinates[1],
                "z": atom.coordinates[2],
                "color": atom.color,
                "display_radius": atom.display_radius,
            }
            for atom in molecule.atoms
        ],
        "bonds": [
            {
                "left": bond.left,
                "right": bond.right,
                "kind": bond.kind,
            }
            for bond in molecule.bonds
        ],
        "dihedrals": [
            {
                "id": dihedral.id,
                "row_index": dihedral.row_index,
                "atom_indices": list(dihedral.atom_indices),
                "atom_labels": list(dihedral.atom_labels),
                "value_degrees": dihedral.value_degrees,
                "kind": dihedral.kind,
                "links": [list(pair) for pair in dihedral.links],
            }
            for dihedral in molecule.dihedrals
        ],
        "warnings": list(molecule.warnings),
    }


HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>
  <style>
    :root {
      color-scheme: light;
      --bg: #f4f6f8;
      --panel: #ffffff;
      --text: #1d2430;
      --muted: #5d6978;
      --line: #d8dee7;
      --line-strong: #b7c0cc;
      --accent: #d76b00;
      --accent-soft: #fff0df;
      --atom-hover: #0077c8;
      --atom-hover-soft: #e6f3ff;
      --focus: #0089a7;
      --shadow: 0 10px 28px rgba(33, 43, 54, 0.10);
    }

    * {
      box-sizing: border-box;
    }

    html,
    body {
      height: 100%;
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      letter-spacing: 0;
    }

    body {
      min-width: 320px;
    }

    .shell {
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      height: 100vh;
    }

    header {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      min-height: 64px;
      padding: 12px 18px;
      border-bottom: 1px solid var(--line);
      background: rgba(255, 255, 255, 0.92);
    }

    h1,
    h2,
    p {
      margin: 0;
    }

    h1 {
      font-size: clamp(18px, 2vw, 24px);
      line-height: 1.2;
      font-weight: 700;
      overflow-wrap: anywhere;
    }

    .source {
      margin-top: 4px;
      color: var(--muted);
      font-size: 12px;
      overflow-wrap: anywhere;
    }

    .summary {
      display: flex;
      align-items: center;
      gap: 8px;
      flex-wrap: wrap;
      justify-content: flex-end;
      color: var(--muted);
      font-size: 13px;
      white-space: nowrap;
    }

    .summary span {
      padding: 5px 8px;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fafbfc;
    }

    main {
      display: grid;
      grid-template-columns: minmax(0, 1fr) minmax(300px, 380px);
      min-height: 0;
    }

    .viewer {
      position: relative;
      min-width: 0;
      min-height: 0;
      background: #eef2f6;
      border-right: 1px solid var(--line);
    }

    canvas {
      display: block;
      width: 100%;
      height: 100%;
      cursor: grab;
      touch-action: none;
    }

    canvas.is-dragging {
      cursor: grabbing;
    }

    canvas.is-atom-hover {
      cursor: pointer;
    }

    .canvas-actions {
      position: absolute;
      top: 14px;
      left: 14px;
      display: flex;
      gap: 8px;
    }

    button {
      font: inherit;
    }

    .tool-button {
      min-height: 34px;
      padding: 0 12px;
      border: 1px solid var(--line-strong);
      border-radius: 6px;
      background: #ffffff;
      color: var(--text);
      box-shadow: 0 1px 3px rgba(33, 43, 54, 0.10);
      cursor: pointer;
    }

    .tool-button:hover,
    .tool-button:focus-visible {
      border-color: var(--focus);
      outline: none;
    }

    aside {
      display: grid;
      grid-template-rows: auto auto minmax(0, 1fr) auto;
      min-width: 0;
      min-height: 0;
      background: var(--panel);
    }

    .side-header {
      padding: 16px 16px 10px;
      border-bottom: 1px solid var(--line);
    }

    h2 {
      font-size: 17px;
      line-height: 1.3;
    }

    .active-detail {
      min-height: 50px;
      padding: 10px 16px;
      border-bottom: 1px solid var(--line);
      color: var(--muted);
      font-size: 13px;
      line-height: 1.45;
    }

    .active-detail strong {
      color: var(--text);
      font-weight: 700;
    }

    .dihedral-list {
      overflow: auto;
      padding: 8px;
    }

    .dihedral-row {
      display: grid;
      width: 100%;
      grid-template-columns: minmax(0, 1fr) auto;
      gap: 6px 12px;
      align-items: center;
      padding: 10px;
      border: 1px solid transparent;
      border-radius: 6px;
      background: transparent;
      color: var(--text);
      text-align: left;
      cursor: pointer;
    }

    .dihedral-row + .dihedral-row {
      margin-top: 4px;
    }

    .dihedral-row:hover,
    .dihedral-row:focus-visible,
    .dihedral-row.is-active,
    .dihedral-row.is-atom-hover {
      border-color: var(--accent);
      background: var(--accent-soft);
      outline: none;
    }

    .dihedral-row.is-atom-hover {
      border-color: var(--atom-hover);
      background: var(--atom-hover-soft);
    }

    .dihedral-row.is-active {
      border-color: var(--accent);
      background: var(--accent-soft);
    }

    .quartet {
      min-width: 0;
      font-size: 13px;
      font-weight: 700;
      overflow-wrap: anywhere;
    }

    .value {
      color: var(--text);
      font-variant-numeric: tabular-nums;
      font-size: 13px;
      white-space: nowrap;
    }

    .meta {
      grid-column: 1 / -1;
      display: flex;
      gap: 8px;
      flex-wrap: wrap;
      color: var(--muted);
      font-size: 12px;
    }

    .kind {
      color: var(--muted);
    }

    .warnings {
      border-top: 1px solid var(--line);
      padding: 10px 16px;
      color: #8a4c00;
      background: #fff8ef;
      font-size: 12px;
      line-height: 1.4;
    }

    .warnings:empty {
      display: none;
    }

    @media (max-width: 820px) {
      .shell {
        height: auto;
        min-height: 100vh;
      }

      header {
        align-items: flex-start;
        flex-direction: column;
      }

      .summary {
        justify-content: flex-start;
      }

      main {
        grid-template-columns: 1fr;
        grid-template-rows: minmax(360px, 56vh) minmax(320px, auto);
      }

      .viewer {
        border-right: 0;
        border-bottom: 1px solid var(--line);
      }
    }
  </style>
</head>
<body>
  <div class="shell">
    <header>
      <div>
        <h1 id="viewerTitle"></h1>
        <p class="source" id="viewerSource"></p>
      </div>
      <div class="summary">
        <span id="atomCount"></span>
        <span id="bondCount"></span>
        <span id="dihedralCount"></span>
      </div>
    </header>
    <main>
      <section class="viewer" aria-label="Molecule visualization">
        <canvas id="moleculeCanvas"></canvas>
        <div class="canvas-actions">
          <button class="tool-button" type="button" id="resetView">Reset</button>
          <button class="tool-button" type="button" id="clearSelection">Clear</button>
        </div>
      </section>
      <aside>
        <div class="side-header">
          <h2>Dihedral Angles</h2>
        </div>
        <div class="active-detail" id="activeDetail"></div>
        <div class="dihedral-list" id="dihedralList"></div>
        <div class="warnings" id="warnings"></div>
      </aside>
    </main>
  </div>
  <script type="application/json" id="moleculeData">__MOLECULE_DATA__</script>
  <script>
    "use strict";

    const molecule = JSON.parse(document.getElementById("moleculeData").textContent);
    window.__viewerMolecule = molecule;

    const canvas = document.getElementById("moleculeCanvas");
    const ctx = canvas.getContext("2d");
    const list = document.getElementById("dihedralList");
    const activeDetail = document.getElementById("activeDetail");
    const buttonsById = new Map();

    const state = {
      rotX: -0.45,
      rotY: 0.72,
      zoom: 1.0,
      hovered: null,
      hoveredAtom: null,
      selected: null,
      dragging: false,
      lastX: 0,
      lastY: 0
    };

    const atomsByIndex = new Map(molecule.atoms.map((atom) => [atom.index, atom]));
    const bondKeys = new Set(molecule.bonds.map((bond) => linkKey(bond.left, bond.right)));
    const center = moleculeCenter(molecule.atoms);
    const radius = Math.max(1.0, moleculeRadius(molecule.atoms, center));

    document.getElementById("viewerTitle").textContent = molecule.title || "Z-matrix viewer";
    document.getElementById("viewerSource").textContent = molecule.source_name || "";
    document.getElementById("atomCount").textContent = molecule.atoms.length + " atoms";
    document.getElementById("bondCount").textContent = molecule.bonds.length + " bonds";
    document.getElementById("dihedralCount").textContent = molecule.dihedrals.length + " dihedrals";
    document.getElementById("warnings").textContent = molecule.warnings.join(" ");

    buildDihedralList();
    resizeCanvas();
    draw();
    window.__viewerReady = true;

    window.addEventListener("resize", () => {
      resizeCanvas();
      draw();
    });

    document.getElementById("resetView").addEventListener("click", () => {
      state.rotX = -0.45;
      state.rotY = 0.72;
      state.zoom = 1.0;
      draw();
    });

    document.getElementById("clearSelection").addEventListener("click", () => {
      state.selected = null;
      state.hovered = null;
      state.hoveredAtom = null;
      canvas.classList.remove("is-atom-hover");
      updateRows();
      draw();
    });

    document.addEventListener("keydown", (event) => {
      if (event.key === "Escape") {
        state.selected = null;
        state.hovered = null;
        state.hoveredAtom = null;
        canvas.classList.remove("is-atom-hover");
        updateRows();
        draw();
      }
    });

    canvas.addEventListener("pointerdown", (event) => {
      state.dragging = true;
      state.lastX = event.clientX;
      state.lastY = event.clientY;
      state.hoveredAtom = null;
      canvas.classList.add("is-dragging");
      canvas.classList.remove("is-atom-hover");
      canvas.setPointerCapture(event.pointerId);
    });

    canvas.addEventListener("pointermove", (event) => {
      if (!state.dragging) {
        const nextHoveredAtom = hitTestAtom(event);
        if (state.hoveredAtom !== nextHoveredAtom) {
          state.hoveredAtom = nextHoveredAtom;
          canvas.classList.toggle("is-atom-hover", state.hoveredAtom !== null);
          updateRows();
          draw();
        }
        return;
      }
      const dx = event.clientX - state.lastX;
      const dy = event.clientY - state.lastY;
      state.lastX = event.clientX;
      state.lastY = event.clientY;
      state.rotY += dx * 0.010;
      state.rotX += dy * 0.010;
      state.rotX = Math.max(-Math.PI * 0.49, Math.min(Math.PI * 0.49, state.rotX));
      draw();
    });

    canvas.addEventListener("pointerup", (event) => {
      state.dragging = false;
      canvas.classList.remove("is-dragging");
      canvas.releasePointerCapture(event.pointerId);
    });

    canvas.addEventListener("pointercancel", () => {
      state.dragging = false;
      state.hoveredAtom = null;
      canvas.classList.remove("is-dragging");
      canvas.classList.remove("is-atom-hover");
      updateRows();
      draw();
    });

    canvas.addEventListener("pointerleave", () => {
      if (state.dragging) {
        return;
      }
      state.hoveredAtom = null;
      canvas.classList.remove("is-atom-hover");
      updateRows();
      draw();
    });

    canvas.addEventListener("wheel", (event) => {
      event.preventDefault();
      const nextZoom = state.zoom * Math.exp(-event.deltaY * 0.001);
      state.zoom = Math.max(0.35, Math.min(5.0, nextZoom));
      draw();
    }, { passive: false });

    function buildDihedralList() {
      if (molecule.dihedrals.length === 0) {
        const empty = document.createElement("p");
        empty.className = "active-detail";
        empty.textContent = "No dihedral rows";
        list.appendChild(empty);
        activeDetail.textContent = "";
        return;
      }

      for (const dihedral of molecule.dihedrals) {
        const row = document.createElement("button");
        row.type = "button";
        row.className = "dihedral-row";
        row.dataset.id = dihedral.id;

        const quartet = document.createElement("span");
        quartet.className = "quartet";
        quartet.textContent = dihedral.atom_labels.join("-");
        row.appendChild(quartet);

        const value = document.createElement("span");
        value.className = "value";
        value.textContent = formatAngle(dihedral.value_degrees);
        row.appendChild(value);

        const meta = document.createElement("span");
        meta.className = "meta";
        meta.textContent = "row " + dihedral.row_index + " / " + dihedral.kind;
        row.appendChild(meta);

        row.addEventListener("mouseenter", () => {
          state.hovered = dihedral.id;
          updateRows();
          draw();
        });
        row.addEventListener("mouseleave", () => {
          state.hovered = null;
          updateRows();
          draw();
        });
        row.addEventListener("focus", () => {
          state.hovered = dihedral.id;
          updateRows();
          draw();
        });
        row.addEventListener("blur", () => {
          state.hovered = null;
          updateRows();
          draw();
        });
        row.addEventListener("click", () => {
          state.selected = state.selected === dihedral.id ? null : dihedral.id;
          updateRows();
          draw();
        });

        buttonsById.set(dihedral.id, row);
        list.appendChild(row);
      }
      updateRows();
    }

    function activeDihedral() {
      const id = state.hovered || state.selected;
      if (!id) {
        return null;
      }
      return molecule.dihedrals.find((dihedral) => dihedral.id === id) || null;
    }

    function atomHoverDihedrals() {
      if (state.hoveredAtom === null) {
        return [];
      }
      return molecule.dihedrals.filter((dihedral) => dihedral.atom_indices[0] === state.hoveredAtom);
    }

    function updateRows() {
      const active = activeDihedral();
      const atomHoverIds = new Set(atomHoverDihedrals().map((dihedral) => dihedral.id));
      for (const [id, row] of buttonsById) {
        row.classList.toggle("is-active", active !== null && active.id === id);
        row.classList.toggle("is-atom-hover", atomHoverIds.has(id) && !(active !== null && active.id === id));
      }
      if (!active) {
        const atomDihedrals = atomHoverDihedrals();
        if (atomDihedrals.length === 0) {
          activeDetail.textContent = "";
          return;
        }
        const atom = atomsByIndex.get(state.hoveredAtom);
        activeDetail.innerHTML = "";
        const strong = document.createElement("strong");
        strong.textContent = atom ? atom.label : "Atom " + state.hoveredAtom;
        activeDetail.appendChild(strong);
        activeDetail.appendChild(document.createTextNode(
          " starts " + atomDihedrals.length + " dihedral" + (atomDihedrals.length === 1 ? "" : "s")
        ));
        return;
      }
      activeDetail.innerHTML = "";
      const strong = document.createElement("strong");
      strong.textContent = active.atom_labels.join("-");
      activeDetail.appendChild(strong);
      activeDetail.appendChild(document.createTextNode(
        "  " + formatAngle(active.value_degrees) + "  " + active.kind
      ));
    }

    function resizeCanvas() {
      const rect = canvas.getBoundingClientRect();
      const dpr = window.devicePixelRatio || 1;
      canvas.width = Math.max(1, Math.floor(rect.width * dpr));
      canvas.height = Math.max(1, Math.floor(rect.height * dpr));
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    function draw() {
      updateRows();
      const width = canvas.clientWidth;
      const height = canvas.clientHeight;
      ctx.clearRect(0, 0, width, height);
      ctx.fillStyle = "#eef2f6";
      ctx.fillRect(0, 0, width, height);

      const active = activeDihedral();
      const atomHover = atomHoverDihedrals();
      const activeAtoms = new Set(active ? active.atom_indices : []);
      const activeLinks = new Set(active ? active.links.map((link) => linkKey(link[0], link[1])) : []);
      const atomHoverAtoms = new Set(atomHover.flatMap((dihedral) => dihedral.atom_indices));
      const atomHoverLinks = new Set(atomHover.flatMap((dihedral) => dihedral.links.map((link) => linkKey(link[0], link[1]))));
      const projected = new Map();
      for (const atom of molecule.atoms) {
        projected.set(atom.index, projectAtom(atom, width, height));
      }

      const bonds = molecule.bonds.map((bond) => {
        const left = projected.get(bond.left);
        const right = projected.get(bond.right);
        return {
          bond,
          left,
          right,
          z: (left.z + right.z) * 0.5
        };
      }).sort((a, b) => a.z - b.z);

      for (const item of bonds) {
        const key = linkKey(item.bond.left, item.bond.right);
        drawBond(item.left, item.right, {
          active: activeLinks.has(key),
          atomHover: atomHoverLinks.has(key) && !activeLinks.has(key),
          inferred: item.bond.kind === "inferred"
        });
      }

      if (atomHover.length > 0) {
        for (const dihedral of atomHover) {
          for (const link of dihedral.links) {
            const key = linkKey(link[0], link[1]);
            if (activeLinks.has(key)) {
              continue;
            }
            const left = projected.get(link[0]);
            const right = projected.get(link[1]);
            drawReferenceLink(left, right, !bondKeys.has(key), "#0077c8", 2.5);
          }
        }
      }

      if (active) {
        for (const link of active.links) {
          const left = projected.get(link[0]);
          const right = projected.get(link[1]);
          drawReferenceLink(left, right, !bondKeys.has(linkKey(link[0], link[1])), "#d76b00", 3);
        }
      }

      const atoms = molecule.atoms.map((atom) => ({
        atom,
        point: projected.get(atom.index)
      })).sort((a, b) => a.point.z - b.point.z);
      for (const item of atoms) {
        drawAtom(
          item.atom,
          item.point,
          activeAtoms.has(item.atom.index),
          atomHoverAtoms.has(item.atom.index) && !activeAtoms.has(item.atom.index)
        );
      }
    }

    function drawBond(left, right, options) {
      ctx.save();
      ctx.beginPath();
      ctx.moveTo(left.x, left.y);
      ctx.lineTo(right.x, right.y);
      ctx.lineCap = "round";
      ctx.lineWidth = options.active ? 8 : options.atomHover ? 7 : 5;
      ctx.strokeStyle = options.active ? "#d76b00" : options.atomHover ? "#0077c8" : options.inferred ? "#b9c3cf" : "#8995a3";
      if (options.inferred && !options.active && !options.atomHover) {
        ctx.setLineDash([6, 6]);
      }
      ctx.stroke();
      ctx.restore();
    }

    function drawReferenceLink(left, right, dashed, color, width) {
      ctx.save();
      ctx.beginPath();
      ctx.moveTo(left.x, left.y);
      ctx.lineTo(right.x, right.y);
      ctx.lineCap = "round";
      ctx.lineWidth = width;
      ctx.strokeStyle = color;
      if (dashed) {
        ctx.setLineDash([7, 5]);
      }
      ctx.stroke();
      ctx.restore();
    }

    function drawAtom(atom, point, highlighted, atomHoverHighlighted) {
      const scale = projectionScale(canvas.clientWidth, canvas.clientHeight);
      const radiusPx = Math.max(5, atom.display_radius * scale);
      if (highlighted || atomHoverHighlighted) {
        ctx.save();
        ctx.beginPath();
        ctx.arc(point.x, point.y, radiusPx + 6, 0, Math.PI * 2);
        ctx.strokeStyle = highlighted ? "#d76b00" : "#0077c8";
        ctx.lineWidth = 4;
        ctx.stroke();
        ctx.restore();
      }

      const gradient = ctx.createRadialGradient(
        point.x - radiusPx * 0.35,
        point.y - radiusPx * 0.45,
        Math.max(1, radiusPx * 0.08),
        point.x,
        point.y,
        radiusPx
      );
      gradient.addColorStop(0, "#ffffff");
      gradient.addColorStop(0.25, atom.color);
      gradient.addColorStop(1, adjustHex(atom.color, -38));

      ctx.save();
      ctx.beginPath();
      ctx.arc(point.x, point.y, radiusPx, 0, Math.PI * 2);
      ctx.fillStyle = gradient;
      ctx.fill();
      ctx.lineWidth = 1.2;
      ctx.strokeStyle = "#26303a";
      ctx.stroke();
      ctx.restore();

      drawLabel(atom.label, point.x + radiusPx + 4, point.y - radiusPx - 2, highlighted, atomHoverHighlighted);
    }

    function drawLabel(text, x, y, highlighted, atomHoverHighlighted) {
      ctx.save();
      ctx.font = (highlighted || atomHoverHighlighted ? "700 " : "600 ") + "12px -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif";
      ctx.lineWidth = 4;
      ctx.strokeStyle = "rgba(255,255,255,0.90)";
      ctx.strokeText(text, x, y);
      ctx.fillStyle = highlighted ? "#a34400" : atomHoverHighlighted ? "#005b99" : "#1d2430";
      ctx.fillText(text, x, y);
      ctx.restore();
    }

    function hitTestAtom(event) {
      const rect = canvas.getBoundingClientRect();
      const x = event.clientX - rect.left;
      const y = event.clientY - rect.top;
      const width = canvas.clientWidth;
      const height = canvas.clientHeight;
      const scale = projectionScale(width, height);
      let best = null;
      let bestDistance = Infinity;

      for (const atom of molecule.atoms) {
        const point = projectAtom(atom, width, height);
        const radiusPx = Math.max(5, atom.display_radius * scale) + 8;
        const dx = x - point.x;
        const dy = y - point.y;
        const distance = Math.sqrt(dx * dx + dy * dy);
        if (distance <= radiusPx && distance < bestDistance) {
          best = atom.index;
          bestDistance = distance;
        }
      }

      return best;
    }

    function projectAtom(atom, width, height) {
      const x0 = atom.x - center.x;
      const y0 = atom.y - center.y;
      const z0 = atom.z - center.z;
      const cy = Math.cos(state.rotY);
      const sy = Math.sin(state.rotY);
      const cx = Math.cos(state.rotX);
      const sx = Math.sin(state.rotX);

      const x1 = x0 * cy + z0 * sy;
      const z1 = -x0 * sy + z0 * cy;
      const y1 = y0 * cx - z1 * sx;
      const z2 = y0 * sx + z1 * cx;
      const scale = projectionScale(width, height);
      return {
        x: width * 0.5 + x1 * scale,
        y: height * 0.5 - y1 * scale,
        z: z2
      };
    }

    function projectionScale(width, height) {
      return Math.min(width, height) / (radius * 2.7) * state.zoom;
    }

    function moleculeCenter(atoms) {
      if (atoms.length === 0) {
        return { x: 0, y: 0, z: 0 };
      }
      const total = atoms.reduce((acc, atom) => ({
        x: acc.x + atom.x,
        y: acc.y + atom.y,
        z: acc.z + atom.z
      }), { x: 0, y: 0, z: 0 });
      return {
        x: total.x / atoms.length,
        y: total.y / atoms.length,
        z: total.z / atoms.length
      };
    }

    function moleculeRadius(atoms, centerPoint) {
      let maxRadius = 0;
      for (const atom of atoms) {
        const dx = atom.x - centerPoint.x;
        const dy = atom.y - centerPoint.y;
        const dz = atom.z - centerPoint.z;
        maxRadius = Math.max(maxRadius, Math.sqrt(dx * dx + dy * dy + dz * dz));
      }
      return maxRadius;
    }

    function linkKey(left, right) {
      return left < right ? left + "-" + right : right + "-" + left;
    }

    function formatAngle(value) {
      return Number(value).toFixed(2) + " deg";
    }

    function adjustHex(hex, amount) {
      const clean = hex.replace("#", "");
      if (clean.length !== 6) {
        return hex;
      }
      const parts = [0, 2, 4].map((index) => {
        const next = Math.max(0, Math.min(255, parseInt(clean.slice(index, index + 2), 16) + amount));
        return next.toString(16).padStart(2, "0");
      });
      return "#" + parts.join("");
    }
  </script>
</body>
</html>
"""
