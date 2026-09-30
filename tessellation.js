/* Draws a faint random spatial structure behind the page header.
   Five models are shown in turn: each visit shows the next one. On the
   Home page, "Resample" draws a new sample of the current model and
   "Next model" switches to the next one.

     1. Poisson–Voronoi tessellation
     2. Boolean model (union of random discs around Poisson germs)
     3. beta-Voronoi tessellation           (Laguerre, beta = 1)
     4. beta'-Voronoi tessellation          (Laguerre, beta = 5)
     5. Gaussian-Voronoi tessellation       (Laguerre)

   The three Laguerre tessellations follow Gusakova, Kabluchko and Thäle:
   a Poisson process of marked points (x, h) in R^2 x R, where the cell of
   (x, h) is the set of y minimising the power |y - x|^2 + h, with height
   intensity  h^beta dh on (0, inf)          for the beta model,
              (-h)^(-beta) dh on (-inf, 0)   for the beta' model,
              e^h dh on R                    for the Gaussian model.
   Heights are truncated where they cannot produce a cell (with overwhelming
   probability), which keeps the simulation fast.

   You never need to edit this file. */
(function () {
  "use strict";

  var STORE_KEY = "tessModel";

  // ---------- random numbers ----------
  var rand = Math.random;
  function poisson(mean) {                 // exact, split into chunks for large means
    var total = 0;
    while (mean > 0) {
      var m = Math.min(mean, 500), L = Math.exp(-m), k = 0, p = 1;
      do { k++; p *= rand(); } while (p > L);
      total += k - 1; mean -= m;
    }
    return total;
  }

  // ---------- geometry ----------
  // Clip convex polygon to { y : 2 y.(b - a) <= |b|^2 - |a|^2 + wb - wa },
  // i.e. the points where the power w.r.t. a is at most the power w.r.t. b.
  function clip(poly, ax, ay, wa, bx, by, wb) {
    var nx = bx - ax, ny = by - ay;
    var c = (bx * bx + by * by - ax * ax - ay * ay + wb - wa) / 2;
    var out = [];
    for (var i = 0, n = poly.length; i < n; i++) {
      var p = poly[i], q = poly[(i + 1) % n];
      var fp = nx * p[0] + ny * p[1] - c, fq = nx * q[0] + ny * q[1] - c;
      if (fp <= 0) out.push(p);
      if ((fp < 0 && fq > 0) || (fp > 0 && fq < 0)) {
        var t = fp / (fp - fq);
        out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
      }
    }
    return out;
  }

  // Laguerre (power) cells of weighted points.  Each point i has position
  // (x[i], y[i]), weight w[i] and reach r[i]: its cell lies inside the disc
  // of radius r[i] around it, and only points with |a - b| < r_a + r_b can
  // share a boundary.  Returns an array of { i, poly }.
  function laguerreCells(x, y, w, r, box) {
    var n = x.length, cells = [];
    for (var a = 0; a < n; a++) {
      if (!(r[a] > 0)) continue;
      var x0 = Math.max(box[0], x[a] - r[a]), x1 = Math.min(box[2], x[a] + r[a]);
      var y0 = Math.max(box[1], y[a] - r[a]), y1 = Math.min(box[3], y[a] + r[a]);
      if (x0 >= x1 || y0 >= y1) continue;
      var poly = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
      for (var b = 0; b < n && poly.length; b++) {
        if (b === a || !(r[b] > 0)) continue;
        var dx = x[b] - x[a], dy = y[b] - y[a], R = r[a] + r[b];
        if (dx * dx + dy * dy >= R * R) continue;
        poly = clip(poly, x[a], y[a], w[a], x[b], y[b], w[b]);
      }
      if (poly.length >= 3) cells.push({ i: a, poly: poly });
    }
    return cells;
  }

  // ---------- the five models ----------
  // All lengths in CSS pixels.  "unit" is the model's length scale in pixels.
  function samplePoints(W, H, m, perUnitArea, unit, height) {
    var x0 = -m, y0 = -m, x1 = W + m, y1 = H + m;
    var areaUnits = (x1 - x0) * (y1 - y0) / (unit * unit);
    var n = poisson(perUnitArea * areaUnits);
    var X = [], Y = [], Hh = [];
    for (var i = 0; i < n; i++) {
      X.push(x0 + rand() * (x1 - x0));
      Y.push(y0 + rand() * (y1 - y0));
      Hh.push(height());                  // height in units^2
    }
    return { x: X, y: Y, h: Hh, box: [x0, y0, x1, y1] };
  }

  // Laguerre model: heights h (units^2), power |y-x|^2 + h, and an upper
  // bound M for the minimal power (holds everywhere with overwhelming probability).
  function laguerreModel(W, H, perUnitArea, unit, height, M) {
    var P = samplePoints(W, H, 3 * unit, perUnitArea, unit, height);
    var u2 = unit * unit, n = P.x.length, w = [], r = [];
    for (var i = 0; i < n; i++) {
      w.push(P.h[i] * u2);
      var rr = (M - P.h[i]) * u2;
      r.push(rr > 0 ? Math.sqrt(rr) : 0);
    }
    return { kind: "cells", x: P.x, y: P.y, cells: laguerreCells(P.x, P.y, w, r, P.box) };
  }

  var MODELS = [
    {
      name: "Poisson–Voronoi tessellation",
      build: function (W, H) {
        // Ordinary Voronoi = Laguerre with all weights 0; reach large enough to be exact here.
        var P = samplePoints(W, H, 60, 1, Math.sqrt(5200), function () { return 0; });
        var n = P.x.length, w = [], r = [];
        for (var i = 0; i < n; i++) { w.push(0); r.push(1e6); }
        return { kind: "cells", x: P.x, y: P.y, cells: laguerreCells(P.x, P.y, w, r, P.box) };
      }
    },
    {
      name: "Boolean model",
      build: function (W, H) {
        // Poisson germs with i.i.d. radii, uniform on [rmin, rmax] (pixels).
        // Intensity chosen for a covered fraction 1 - exp(-lambda * pi * E[R^2]) of about 0.55.
        var m = 50, rmin = 14, rmax = 44, coverage = 0.55;
        var ER2 = (Math.pow(rmax, 3) - Math.pow(rmin, 3)) / (3 * (rmax - rmin));
        var lambda = -Math.log(1 - coverage) / (Math.PI * ER2);
        var x0 = -m, y0 = -m, x1 = W + m, y1 = H + m;
        var n = poisson(lambda * (x1 - x0) * (y1 - y0)), discs = [];
        for (var i = 0; i < n; i++) {
          discs.push([x0 + rand() * (x1 - x0), y0 + rand() * (y1 - y0), rmin + rand() * (rmax - rmin)]);
        }
        return { kind: "discs", discs: discs };
      }
    },
    {
      name: "β-Voronoi tessellation (β = 1)",
      build: function (W, H) {
        // intensity h^1 dh on [0, T]: T = 3 suffices (P(min power > 3) = exp(-9π/2)).
        var T = 3;
        return laguerreModel(W, H, T * T / 2, 68, function () { return T * Math.sqrt(rand()); }, T);
      }
    },
    {
      name: "β′-Voronoi tessellation (β = 5)",
      build: function (W, H) {
        // intensity (-h)^(-5) dh on (-inf, -eps]: count per unit area eps^(-4)/4;
        // P(min power > -eps) = exp(-pi eps^(-3)/12), negligible for eps = 0.3.
        var beta = 5, eps = 0.3;
        return laguerreModel(W, H, Math.pow(eps, 1 - beta) / (beta - 1), 114,
          function () { return -eps * Math.pow(1 - rand(), -1 / (beta - 1)); }, -eps);
      }
    },
    {
      name: "Gaussian–Voronoi tessellation",
      build: function (W, H) {
        // intensity e^h dh on (-inf, Hmax]: count per unit area e^Hmax.
        var Hmax = 2;
        return laguerreModel(W, H, Math.exp(Hmax), 51,
          function () { return Hmax + Math.log(1 - rand()); }, Hmax);
      }
    }
  ];

  // ---------- drawing ----------
  function inView(poly, W, H) {
    var a = Infinity, b = Infinity, c = -Infinity, d = -Infinity;
    for (var k = 0; k < poly.length; k++) {
      a = Math.min(a, poly[k][0]); c = Math.max(c, poly[k][0]);
      b = Math.min(b, poly[k][1]); d = Math.max(d, poly[k][1]);
    }
    return c > 0 && a < W && d > 0 && b < H;
  }

  function draw(canvas, idx) {
    var box = canvas.parentElement.getBoundingClientRect();
    var W = Math.max(1, Math.round(box.width)), H = Math.max(1, Math.round(box.height));
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = W * dpr; canvas.height = H * dpr;
    canvas.style.width = W + "px"; canvas.style.height = H + "px";
    var ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);

    var colour = getComputedStyle(document.documentElement).getPropertyValue("--accent").trim() || "#1d4f91";
    ctx.strokeStyle = colour; ctx.fillStyle = colour;
    ctx.lineWidth = 0.9; ctx.lineJoin = "round";

    var model = MODELS[idx], S = model.build(W, H), count = 0;

    if (S.kind === "cells") {
      for (var c = 0; c < S.cells.length; c++) {
        var poly = S.cells[c].poly, i = S.cells[c].i;
        if (inView(poly, W, H)) count++;
        ctx.globalAlpha = 0.34;
        ctx.beginPath();
        ctx.moveTo(poly[0][0], poly[0][1]);
        for (var k = 1; k < poly.length; k++) ctx.lineTo(poly[k][0], poly[k][1]);
        ctx.closePath(); ctx.stroke();
        ctx.globalAlpha = 0.55;              // nucleus (may lie outside its Laguerre cell)
        ctx.beginPath(); ctx.arc(S.x[i], S.y[i], 1.6, 0, 2 * Math.PI); ctx.fill();
      }
    } else {
      var D = S.discs;
      // union, filled once so overlaps do not look darker
      var off = document.createElement("canvas");
      off.width = canvas.width; off.height = canvas.height;
      var o = off.getContext("2d");
      o.setTransform(dpr, 0, 0, dpr, 0, 0); o.fillStyle = colour;
      for (var j = 0; j < D.length; j++) { o.beginPath(); o.arc(D[j][0], D[j][1], D[j][2], 0, 2 * Math.PI); o.fill(); }
      ctx.save(); ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.globalAlpha = 0.1; ctx.drawImage(off, 0, 0); ctx.restore();
      // boundary of the union: arcs of each circle not covered by another disc
      ctx.globalAlpha = 0.4;
      for (j = 0; j < D.length; j++) {
        var steps = 96, prev = null;
        for (var s = 0; s <= steps; s++) {
          var th = 2 * Math.PI * s / steps;
          var px = D[j][0] + D[j][2] * Math.cos(th), py = D[j][1] + D[j][2] * Math.sin(th), covered = false;
          for (var q = 0; q < D.length && !covered; q++) {
            if (q === j) continue;
            var ex = px - D[q][0], ey = py - D[q][1];
            covered = ex * ex + ey * ey < D[q][2] * D[q][2];
          }
          var pt = covered ? null : [px, py];
          if (prev && pt) { ctx.beginPath(); ctx.moveTo(prev[0], prev[1]); ctx.lineTo(pt[0], pt[1]); ctx.stroke(); }
          prev = pt;
        }
        if (D[j][0] > 0 && D[j][0] < W && D[j][1] > 0 && D[j][1] < H) count++;
        ctx.globalAlpha = 0.55;
        ctx.beginPath(); ctx.arc(D[j][0], D[j][1], 1.4, 0, 2 * Math.PI); ctx.fill();
        ctx.globalAlpha = 0.4;
      }
    }
    ctx.globalAlpha = 1;

    var label = document.querySelector("[data-tess-label]");
    if (label) label.textContent = model.name;
  }

  // ---------- setup ----------
  function loadIndex() {
    var i = -1;
    try { i = parseInt(window.localStorage.getItem(STORE_KEY), 10); } catch (e) {}
    if (!(i >= 0)) i = Math.floor(rand() * MODELS.length) - 1;
    i = (i + 1) % MODELS.length;
    try { window.localStorage.setItem(STORE_KEY, String(i)); } catch (e) {}
    return i;
  }

  function init() {
    var canvas = document.querySelector(".tess");
    if (!canvas || !canvas.getContext) return;
    var idx = loadIndex();

    // Home page: replace the caption button by a caption and two buttons.
    var old = document.querySelector("[data-tess-resample]");
    if (old) {
      if (!document.getElementById("tess-style")) {
        var st = document.createElement("style");
        st.id = "tess-style";
        st.textContent =
          ".tess-note .tess-bar{display:flex;flex-wrap:wrap;justify-content:flex-end;align-items:baseline;gap:.15rem .1rem;" +
          "background:var(--paper);background:color-mix(in srgb,var(--paper) 82%,transparent);border-radius:2px;padding:.15rem .5rem;margin-right:-.5rem}" +
          ".tess-note .tess-cap{font:italic var(--step--1)/1.4 var(--serif);color:var(--muted)}" +
          ".tess-note .tess-bar button{font:var(--step--1)/1.4 var(--serif);font-style:normal;color:var(--accent);background:none;" +
          "border:0;padding:0 .3rem;margin:0;cursor:pointer}" +
          ".tess-note .tess-bar button:hover{text-decoration:underline;text-underline-offset:.18em}" +
          ".tess-note .tess-bar button:focus-visible{outline:2px solid var(--accent);outline-offset:2px;border-radius:2px}" +
          ".tess-note .tess-sep{color:var(--muted);padding:0 .15rem}" +
          ".tess-note .tess-acts{white-space:nowrap;display:inline-flex;align-items:baseline}";
        document.head.appendChild(st);
      }
      var bar = document.createElement("div");
      bar.className = "tess-bar";
      bar.innerHTML = '<span class="tess-cap" data-tess-label>…</span>' +
        '<span class="tess-acts"><span class="tess-sep" aria-hidden="true">·</span>' +
        '<button type="button" data-act="resample" title="Draw a new sample of this model">Resample</button>' +
        '<span class="tess-sep" aria-hidden="true">·</span>' +
        '<button type="button" data-act="next" title="Show the next model">Next model</button></span>';
      old.parentNode.replaceChild(bar, old);
      bar.addEventListener("click", function (e) {
        var b = e.target.closest ? e.target.closest("button") : null;
        if (!b) return;
        if (b.getAttribute("data-act") === "next") {
          idx = (idx + 1) % MODELS.length;
          try { window.localStorage.setItem(STORE_KEY, String(idx)); } catch (err) {}
        }
        draw(canvas, idx);
      });
    }
    draw(canvas, idx);

    var timer, lastW = canvas.parentElement.clientWidth;
    window.addEventListener("resize", function () {
      clearTimeout(timer);
      timer = setTimeout(function () {
        var wNow = canvas.parentElement.clientWidth;
        if (wNow !== lastW) { lastW = wNow; draw(canvas, idx); }
      }, 200);
    });
    var redraw = function () { draw(canvas, idx); };
    if (window.matchMedia) {
      var mq = window.matchMedia("(prefers-color-scheme: dark)");
      if (mq.addEventListener) mq.addEventListener("change", redraw);
    }
    new MutationObserver(redraw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
