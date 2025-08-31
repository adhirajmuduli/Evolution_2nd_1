// static/script.js
document.addEventListener("DOMContentLoaded", () => {
    // Add this new code at the top of your existing DOMContentLoaded event listener
    const themeToggle = document.getElementById("theme-toggle");
    const currentTheme = localStorage.getItem("theme");

    if (currentTheme === "dark") {
    document.body.classList.add("dark-mode");
    themeToggle.innerHTML = "🌙";
    } else {
    themeToggle.innerHTML = "☀️";
    }

    themeToggle.addEventListener("click", () => {
    document.body.classList.toggle("dark-mode");
    let theme = "light";
    if (document.body.classList.contains("dark-mode")) {
        theme = "dark";
        themeToggle.innerHTML = "🌙";
    } else {
        themeToggle.innerHTML = "☀️";
    }
    localStorage.setItem("theme", theme);
    });
    const numAlleles = document.getElementById("numAlleles");
    const pInput = document.getElementById("p");
    const qInput = document.getElementById("q");
    const rInput = document.getElementById("r");
    const rField = document.getElementById("rField");
  
    const computeBtn = document.getElementById("computeBtn");
    const resultsPre = document.getElementById("results");
  
    const derivationDiv = document.getElementById("derivation");
  
    const genotypeSelect = document.getElementById("genotypeSelect");
    const resolutionInput = document.getElementById("resolution");
    const plotBtn = document.getElementById("plotBtn");
    const plotDiv = document.getElementById("plot");
    const simplexBtn = document.getElementById("simplexBtn");
    const simplexDiv = document.getElementById("simplex");
    // Toggle r input visibility for 2 vs 3 alleles
    function setAlleleUI() {
      if (numAlleles.value === "2") {
        rField.style.display = "none";
      } else {
        rField.style.display = "";
      }
      loadDerivation(); // update algebraic derivation
    }
    numAlleles.addEventListener("change", setAlleleUI);
    setAlleleUI();

    function markInvalid(input, valid) {
        if (valid) {
          input.style.borderColor = "";
        } else {
          input.style.borderColor = "red";
        }
      }      
  
    // Fetch and render detailed symbolic derivation using MathJax
    async function loadDerivation() {
        const n = numAlleles.value;
        try {
        const res = await fetch(`/derive?num_alleles=${n}&detailed=true`);
        const data = await res.json();
        if (data.error) {
            derivationDiv.innerHTML = `<pre style="color:red">${data.error}</pre>`;
            return;
        }
        
        const steps = data.derivation.derivation_steps;
        let html = "";
        
        for (const step of steps) {
            html += `<div class="derivation-step">`;
            html += `<h4>Step ${step.step}: ${step.title}</h4>`;
            html += `<p>${step.content}</p>`;
            
            if (step.constraint) {
            html += `<p style="font-style: italic; color: #666;">\\(${step.constraint}\\)</p>`;
            }
            
            if (step.equations) {
            html += `<ul>`;
            for (const eq of step.equations) {
                html += `<li>\\(${eq}\\)</li>`;
            }
            html += `</ul>`;
            }
            
            if (step.calculation) {
            html += `<p>\\(${step.calculation}\\)</p>`;
            }
            
            if (step.calculations) {
            for (const calc of step.calculations) {
                html += `<div class="allele-calc">`;
                html += `<h5>For allele ${calc.allele}:</h5>`;
                html += `<p>\\(${calc.formula}\\)</p>`;
                html += `<p><strong>Substituting:</strong></p><ul>`;
                for (const sub of calc.substitution) {
                html += `<li>\\(${sub}\\)</li>`;
                }
                html += `</ul>`;
                html += `<p>\\(${calc.expanded}\\)</p>`;
                html += `<p>\\(${calc.simplified}\\)</p>`;
                html += `<p><strong>∴ \\(${calc.conclusion}\\)</strong></p>`;
                html += `</div>`;
            }
            }
            
            if (step.conclusion) {
            html += `<ul>`;
            for (const conc of step.conclusion) {
                html += `<li>\\(${conc}\\)</li>`;
            }
            html += `</ul>`;
            }
            
            if (step.final) {
            html += `<p style="font-weight: bold; color: #1158c7;">${step.final}</p>`;
            }
            
            html += `</div><hr style="margin: 20px 0; border: none; border-top: 1px solid #eee;">`;
        }
        
        derivationDiv.innerHTML = html;
    
        if (window.MathJax && MathJax.typesetPromise) {
            MathJax.typesetPromise([derivationDiv]).catch((err) => console.error(err));
        }
        } catch (err) {
        console.error(err);
        derivationDiv.innerText = "Error loading derivation";
        }
    }
  
    // Normalize helper: normalize p,q,(r) so sum -> 1
    function normalizeInputs(p, q, r, want3) {
      if (want3) {
        let s = p + q + r;
        if (s <= 0) throw new Error("At least one allele freq must be positive");
        if (Math.abs(s - 1) > 1e-8) {
          return [p / s, q / s, r / s, true];
        }
        return [p, q, r, false];
      } else {
        let s = p + q;
        if (s <= 0) throw new Error("At least one allele freq must be positive");
        if (Math.abs(s - 1) > 1e-8) {
          return [p / s, q / s, 0.0, true];
        }
        return [p, q, 0.0, false];
      }
    }
  
    // Compute numeric genotype frequencies and show conservation steps
    // Compute numeric genotype frequencies and show conservation steps
    computeBtn.addEventListener("click", async () => {
        try {
          let p = parseFloat(pInput.value) || 0;
          let q = parseFloat(qInput.value) || 0;
          let r = parseFloat(rInput.value) || 0;
          const want3 = numAlleles.value === "3";
  
          // quick client-side check
          const total = p + q + (want3 ? r : 0);
          if (Math.abs(total - 1) > 0.01) {
          alert(`Warning: p+q${want3 ? "+r" : ""} = ${total.toFixed(3)} (should be close to 1). Values will be normalized automatically.`);
          }
    
          const [pn, qn, rn, normalized] = normalizeInputs(p, q, r, want3);
          // update inputs to normalized values
          pInput.value = pn.toFixed(6);
          qInput.value = qn.toFixed(6);
          markInvalid(pInput, p >= 0 && p <= 1);
          markInvalid(qInput, q >= 0 && q <= 1);
          if (want3) markInvalid(rInput, r >= 0 && r <= 1);
          if (want3) rInput.value = rn.toFixed(6);
    
          // call backend
          const params = new URLSearchParams({ p: pn, q: qn, r: rn });
          const resp = await fetch(`/api/genotypes?${params.toString()}`);
          const data = await resp.json();
  
          if (data.error) {
          resultsPre.textContent = `Error: ${data.error}`;
          return;
          }
    
          // Nicely format numeric results and steps
          const lines = [];
          lines.push("Input allele freqs:");
          lines.push(`  p = ${data.input.p.toFixed(6)}, q = ${data.input.q.toFixed(6)}, r = ${data.input.r.toFixed(6)}`);
          lines.push("");
          lines.push("Genotype frequencies:");
          for (const [g, v] of Object.entries(data.verbose_genotypes)) {
            lines.push(`  ${g}: ${v}`);
          }
          lines.push("");
          lines.push("Allele frequency reconstruction (next generation):");
          lines.push(`  p' = ${data.conservation.p_next.value.toFixed(8)}  (${data.conservation.p_next.formula} -> ${data.conservation.p_next.substitution})`);
          lines.push(`  q' = ${data.conservation.q_next.value.toFixed(8)}  (${data.conservation.q_next.formula} -> ${data.conservation.q_next.substitution})`);
          lines.push(`  r' = ${data.conservation.r_next.value.toFixed(8)}  (${data.conservation.r_next.formula} -> ${data.conservation.r_next.substitution})`);
          resultsPre.textContent = lines.join("\n");
        } catch (err) {
          resultsPre.textContent = `Error: ${err.message}`;
        }
    });
  
    // Plot surface for a selected genotype by consuming /api/grid
    plotBtn.addEventListener("click", async () => {
        const selected = Array.from(genotypeSelect.selectedOptions).map(o => o.value);
        const n = parseInt(resolutionInput.value) || 60;
    
        plotDiv.innerHTML = "<p style='color:gray'>Loading… please wait</p>";
    
        try {
        const traces = [];
    
        for (const genotype of selected) {
            const resp = await fetch(`/api/grid?genotype=${genotype}&n=${n}`);
            const data = await resp.json();
            if (data.error) continue;
    
            // Find the min/max frequencies in the z_matrix
            const z_flattened = data.z_matrix.flat().filter(z => z !== null);
            const cmin = Math.min(...z_flattened);
            const cmax = Math.max(...z_flattened);

            traces.push({
                x: data.p_vals,
                y: data.q_vals,
                z: data.z_matrix,
                type: 'surface',
                name: genotype,
                showscale: traces.length === 0,  // show only one colorbar
                colorscale: "turbo",
                cmin: cmin, // <-- Updated line
                cmax: cmax, // <-- Updated line
                colorbar: {
                    title: 'Frequency',
                    titleside: 'right'
                },
                contours: { // contour projection
                    z: {
                        show: true,
                        usecolormap: true,
                        highlightcolor: "#42f462",
                        project: { z: true }
                    },
                    x: { show: true, color: "black", width: 1 },
                    y: { show: true, color: "black", width: 1 },
                }
            });
        }
    
        if (selected.length === 0) {
            plotDiv.textContent = "Select at least one genotype to plot.";
            return;
        }
    
        const layout = {
            title: {
            text: `Genotype surfaces: ${selected.join(", ")}`,
            font: { size: 20 }
            },
            scene: {
            xaxis: { title: 'p (allele A freq)', range: [0, 1] },
            yaxis: { title: 'q (allele B freq)', range: [0, 1] },
            zaxis: { title: 'Genotype Frequency', range: [0, 1] },
            camera: {
                eye: { x: 1.6, y: 1.6, z: 0.9 }   // better angle
            },
            aspectmode: "cube"                  // equal scaling
            },
            autosize: true,
            margin: { l: 0, r: 0, b: 0, t: 40 }
        };
    
        plotDiv.innerHTML = "";
        Plotly.newPlot(plotDiv, traces, layout, {responsive: true});
        } catch (err) {
        plotDiv.innerHTML = `<p style="color:red">Error plotting: ${err.message}</p>`;
        console.error(err);
        }
    });
  

// Plot 2D simplex (ternary plot) for one genotype
simplexBtn.addEventListener("click", async () => {
    const selected = Array.from(genotypeSelect.selectedOptions).map(o => o.value);
    if (selected.length !== 1) {
    simplexDiv.textContent = "Please select exactly one genotype for simplex plot.";
    return;
    }
    const genotype = selected[0];

    const n = parseInt(resolutionInput.value) || 60;
  
    try {
      const resp = await fetch(`/api/grid?genotype=${genotype}&n=${n}`);
      const data = await resp.json();
      if (data.error) {
        simplexDiv.textContent = `Error: ${data.error}`;
        return;
      }
  
      const pts = data.points;
      if (!pts || pts.length === 0) {
        simplexDiv.textContent = "No grid points returned.";
        return;
      }
  
      // Replace the simplexTrace code with this:
      // Calculate min and max frequencies from the data
      const freqs = pts.map(pt => pt.freq).filter(f => f !== null);
      const cmin = Math.min(...freqs);
      const cmax = Math.max(...freqs);
      console.log("Frequency range:", cmin, cmax);
      console.log("Number of points:", pts.length);

      const simplexTrace = {
        type: 'scatterternary',
        mode: 'markers',
        a: pts.map(pt => pt.p ),
        b: pts.map(pt => pt.q ),
        c: pts.map(pt => pt.r ),
        marker: {
            color: pts.map(pt => pt.freq),
            colorscale: 'turbo',
            cmin: cmin, // <-- Updated line
            cmax: cmax, // <-- Updated line
            size: 6,
            showscale: true,
            colorbar: {
                title: `${genotype} Frequency`,
                titleside: 'right'
            },
            opacity: 1
        },
        text: pts.map(pt =>
            `p: ${pt.p.toFixed(3)}<br>q: ${pt.q.toFixed(3)}<br>r: ${pt.r.toFixed(3)}<br>${genotype}: ${pt.freq.toFixed(4)}`
        ),
        hoverinfo: 'text',
        name: genotype
      };
      
      // Replace the existing simplexLayout code with this:
      const simplexLayout = {
        title: `Genotype frequency simplex (${genotype})`,
        ternary: {
            sum: 1,
            aaxis: {
                title: { 
                    text: 'p', 
                    font: { size: 16, weight: 'bold' }
                },
                min: 0,
                linewidth: 2,
                ticks: 'outside',
                tickangle: 0,
                tickfont: { size: 10 },
                color: '#1f77b4',
                hoverformat: '.0f'
            },
            baxis: {
                title: { 
                    text: 'q', 
                    font: { size: 16, weight: 'bold' }
                },
                min: 0,
                linewidth: 2,
                ticks: 'outside',
                tickangle: 0,
                tickfont: { size: 10 },
                color: '#ff7f0e',
                hoverformat: '.0f'
            },
            caxis: {
                title: { 
                    text: 'r', 
                    font: { size: 16, weight: 'bold' }
                },
                min: 0,
                linewidth: 2,
                ticks: 'outside',
                tickangle: 0,
                tickfont: { size: 10 },
                color: '#2ca02c',
                hoverformat: '.0f'
            }
        },
        showlegend: false,
        annotations: [
            {
                text: '.',
                x: 0.5,
                y: 1.05,
                xref: 'paper',
                yref: 'paper',
                showarrow: false,
                font: { size: 12, color: '#1f77b4' }
            },
            {
                text: '.',
                x: 0.95,
                y: 0.05,
                xref: 'paper',
                yref: 'paper',
                showarrow: false,
                font: { size: 12, color: '#ff7f0e' }
            },
            {
                text: '.',
                x: 0.05,
                y: 0.05,
                xref: 'paper',
                yref: 'paper',
                showarrow: false,
                font: { size: 12, color: '#2ca02c' }
            }
        ],
        margin: { t: 80, b: 60, l: 60, r: 60 }
      };
  
      simplexDiv.innerHTML = "";
      Plotly.newPlot(simplexDiv, [simplexTrace], simplexLayout);
    } catch (err) {
      simplexDiv.textContent = `Error plotting simplex: ${err.message}`;
      console.error(err);
    }
  });  
  
    // initial load
    loadDerivation();
  });
  