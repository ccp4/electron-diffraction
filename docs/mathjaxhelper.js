MathJax.Hub.Config({
  config: ["MMLorHTML.js"],
  jax: ["input/TeX", "output/HTML-CSS", "output/NativeMML"],
  extensions: ["MathMenu.js", "MathZoom.js"],
  TeX: {
    Macros: {
        cc: ["{\\mathcal #1}",1],
        bb: ["{\\mathbf #1}",1],
        dP: "{\\partial}",
        grad: "{\\nabla}"
    }
  }
});
