/* Draws a faint Poisson–Voronoi tessellation behind the page header.
   A new sample is drawn on each visit (and when the caption is clicked).
   You never need to edit this file. */
(function () {
  "use strict";

  var INTENSITY = 1 / 5200;   // expected number of nuclei per square pixel

  function poisson(mean, rand) {           // Knuth's method; fine for mean < 700
    var L = Math.exp(-mean), k = 0, p = 1;
    do { k++; p *= rand(); } while (p > L);
    return k - 1;
  }

  // Clip convex polygon by the half-plane of points closer to a than to b.
  function clip(poly, a, b) {
    var nx = b[0] - a[0], ny = b[1] - a[1];
    var c = (b[0] * b[0] + b[1] * b[1] - a[0] * a[0] - a[1] * a[1]) / 2;
    var out = [];
    for (var i = 0; i < poly.length; i++) {
      var p = poly[i], q = poly[(i + 1) % poly.length];
      var fp = nx * p[0] + ny * p[1] - c, fq = nx * q[0] + ny * q[1] - c;
      if (fp <= 0) out.push(p);
      if ((fp < 0 && fq > 0) || (fp > 0 && fq < 0)) {
        var t = fp / (fp - fq);
        out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
      }
    }
    return out;
  }

  function draw(canvas) {
    var box = canvas.parentElement.getBoundingClientRect();
    var W = Math.max(1, Math.round(box.width)), H = Math.max(1, Math.round(box.height));
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = W * dpr; canvas.height = H * dpr;
    canvas.style.width = W + "px"; canvas.style.height = H + "px";
    var ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);

    var css = getComputedStyle(document.documentElement);
    var colour = css.getPropertyValue("--accent").trim() || "#1d4f91";

    // Sample on a slightly larger window so cells near the edge are not cut short.
    var m = 60, x0 = -m, y0 = -m, x1 = W + m, y1 = H + m;
    var rand = Math.random;
    var n = poisson(INTENSITY * (x1 - x0) * (y1 - y0), rand);
    var pts = [];
    for (var i = 0; i < n; i++) pts.push([x0 + rand() * (x1 - x0), y0 + rand() * (y1 - y0)]);

    ctx.strokeStyle = colour;
    ctx.fillStyle = colour;
    ctx.lineWidth = 0.9;
    ctx.lineJoin = "round";

    for (var a = 0; a < n; a++) {
      var cell = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]];
      for (var b = 0; b < n && cell.length; b++) {
        if (a !== b) cell = clip(cell, pts[a], pts[b]);
      }
      if (cell.length < 3) continue;
      ctx.globalAlpha = 0.34;
      ctx.beginPath();
      ctx.moveTo(cell[0][0], cell[0][1]);
      for (var k = 1; k < cell.length; k++) ctx.lineTo(cell[k][0], cell[k][1]);
      ctx.closePath();
      ctx.stroke();
      ctx.globalAlpha = 0.55;
      ctx.beginPath();
      ctx.arc(pts[a][0], pts[a][1], 1.6, 0, 2 * Math.PI);
      ctx.fill();
    }
    ctx.globalAlpha = 1;

    var count = document.querySelector("[data-tess-count]");
    if (count) {
      var visible = 0;
      for (var v = 0; v < n; v++) if (pts[v][0] >= 0 && pts[v][0] <= W && pts[v][1] >= 0 && pts[v][1] <= H) visible++;
      count.textContent = String(visible);
    }
  }

  function init() {
    var canvas = document.querySelector(".tess");
    if (!canvas || !canvas.getContext) return;
    draw(canvas);

    var timer, lastW = canvas.parentElement.clientWidth;
    window.addEventListener("resize", function () {
      clearTimeout(timer);
      timer = setTimeout(function () {
        var w = canvas.parentElement.clientWidth;
        if (w !== lastW) { lastW = w; draw(canvas); }
      }, 200);
    });

    var redraw = function () { draw(canvas); };
    if (window.matchMedia) {
      var mq = window.matchMedia("(prefers-color-scheme: dark)");
      if (mq.addEventListener) mq.addEventListener("change", redraw);
    }
    new MutationObserver(redraw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });

    var btn = document.querySelector("[data-tess-resample]");
    if (btn) btn.addEventListener("click", redraw);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
