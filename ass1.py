# backend/app.py
from flask import Flask, request, jsonify
import numpy as np
import csv
from flask import Response 
import sympy
import os
import logging
import math
import werkzeug
from sympy import symbols, Eq, solve, expand, simplify, latex
app = Flask(__name__)

# -------------------------------
# Core HW logic (3 alleles, 1 locus)
# -------------------------------
from flask import render_template

@app.route("/")
def index():
    return render_template("ass1.html")

def safe(val):
    if isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
        return None  # or 0.0 depending on what makes sense
    return val

def genotype_frequencies(*alleles):
    """
    Hardy-Weinberg genotype frequencies for scalar allele frequencies.
    """
    if abs(sum(alleles) - 1.0) > 1e-8:
        raise ValueError("Allele frequencies must sum to 1 (within tolerance).")

    freqs = {}
    for i, pi in enumerate(alleles):
        freqs[f"A{i+1}A{i+1}"] = pi**2
        for j, pj in enumerate(alleles[i+1:], i+1):
            freqs[f"A{i+1}A{j+1}"] = 2 * pi * pj
    return freqs

def genotype_grid(n=50, genotype="A1A1"):
    p_vals = np.linspace(0, 1, n)
    q_vals = np.linspace(0, 1, n)
    P, Q = np.meshgrid(p_vals, q_vals)
    R = 1 - P - Q
    
    # Apply constraint mask
    mask = R >= 0
    
    if genotype == "A1A1":
        F = np.where(mask, P**2, np.nan)
    elif genotype == "A2A2":
        F = np.where(mask, Q**2, np.nan)
    elif genotype == "A3A3":
        F = np.where(mask, R**2, np.nan)
    elif genotype in ("A1A2","A2A1"):
        F = np.where(mask, 2*P*Q, np.nan)
    elif genotype in ("A1A3","A3A1"):
        F = np.where(mask, 2*P*R, np.nan)
    elif genotype in ("A2A3","A3A2"):
        F = np.where(mask, 2*Q*R, np.nan)
    else:
        raise ValueError("Invalid genotype")
    
    return {
        "p_vals": p_vals.tolist(),
        "q_vals": q_vals.tolist(), 
        "z_matrix": F.tolist()
    }

def derive_hw_equilibrium_detailed(num_alleles=3):
    """
    Generate a comprehensive step-by-step Hardy-Weinberg derivation
    that shows conservation of allele frequencies across generations.
    """
    import sympy
    from sympy import symbols, latex, simplify, expand
    
    if num_alleles == 2:
        p, q = symbols('p q', real=True, positive=True)
        alleles = [p, q]
        allele_names = ['A', 'B']
        constraint = "p + q = 1"
    elif num_alleles == 3:
        p, q, r = symbols('p q r', real=True, positive=True)
        alleles = [p, q, r]
        allele_names = ['A', 'B', 'C']
        constraint = "p + q + r = 1"
    else:
        raise ValueError("Only 2 or 3 alleles supported")
    
    derivation_steps = []
    
    # Step 1: Initial setup
    derivation_steps.append({
        "step": 1,
        "title": "Initial Setup",
        "content": f"Let allele frequencies be: {', '.join([f'f({name}) = {latex(freq)}' for name, freq in zip(allele_names, alleles)])}",
        "constraint": f"Constraint: {latex(sum(alleles))} = 1"
    })
    
    # Step 2: Hardy-Weinberg genotype frequencies
    genotype_freqs = {}
    hw_equations = []
    
    # Homozygotes
    for i, (name, freq) in enumerate(zip(allele_names, alleles)):
        genotype = f"{name}{name}"
        genotype_freqs[genotype] = freq**2
        hw_equations.append(f"f({genotype}) = {latex(freq)}^2 = {latex(freq**2)}")
    
    # Heterozygotes
    for i in range(len(alleles)):
        for j in range(i+1, len(alleles)):
            name1, name2 = allele_names[i], allele_names[j]
            freq1, freq2 = alleles[i], alleles[j]
            genotype = f"{name1}{name2}"
            genotype_freqs[genotype] = 2 * freq1 * freq2
            hw_equations.append(f"f({genotype}) = 2 \\cdot {latex(freq1)} \\cdot {latex(freq2)} = {latex(2 * freq1 * freq2)}")
    
    derivation_steps.append({
        "step": 2,
        "title": "Hardy-Weinberg Genotype Frequencies",
        "content": "Under random mating, genotype frequencies are:",
        "equations": hw_equations
    })
    
    # Step 3: Verification that frequencies sum to 1
    total_freq = sum(genotype_freqs.values())
    expanded_total = expand(total_freq)
    
    derivation_steps.append({
        "step": 3,
        "title": "Verification: Genotype Frequencies Sum to 1",
        "content": f"Sum of all genotype frequencies:",
        "calculation": f"{latex(total_freq)} = {latex(expanded_total)} = {latex(sum(alleles))}^2 = 1^2 = 1"
    })
    
    # Step 4: Next generation allele frequency calculations (the key conservation proof)
    next_gen_calcs = []
    
    for i, (allele_name, original_freq) in enumerate(zip(allele_names, alleles)):
        # Calculate f(allele)' from genotype frequencies
        reconstruction_terms = []
        reconstruction_latex_terms = []
        
        # Homozygote contribution (full frequency)
        homo_genotype = f"{allele_name}{allele_name}"
        homo_contrib = genotype_freqs[homo_genotype]
        reconstruction_terms.append(homo_contrib)
        reconstruction_latex_terms.append(f"f({homo_genotype})")
        
        # Heterozygote contributions (half frequency each)
        for j, other_name in enumerate(allele_names):
            if i != j:
                # Check both orders (AB and BA are the same)
                genotype1 = f"{allele_name}{other_name}"
                genotype2 = f"{other_name}{allele_name}"
                
                if genotype1 in genotype_freqs:
                    hetero_contrib = genotype_freqs[genotype1] / 2
                    reconstruction_terms.append(hetero_contrib)
                    reconstruction_latex_terms.append(f"\\frac{1}{2}f({genotype1})")
                elif genotype2 in genotype_freqs:
                    hetero_contrib = genotype_freqs[genotype2] / 2
                    reconstruction_terms.append(hetero_contrib)
                    reconstruction_latex_terms.append(f"\\frac{1}{2}f({genotype2})")
        
        # Sum all contributions
        next_freq = sum(reconstruction_terms)
        simplified_next_freq = simplify(next_freq)
        
        # Detailed step-by-step calculation
        substitution_steps = []
        for term, latex_term in zip(reconstruction_terms, reconstruction_latex_terms):
            substitution_steps.append(f"{latex_term} = {latex(term)}")
        
        calculation_detail = {
            "allele": allele_name,
            "formula": f"f({allele_name})' = {' + '.join(reconstruction_latex_terms)}",
            "substitution": substitution_steps,
            "expanded": f"f({allele_name})' = {latex(next_freq)}",
            "simplified": f"f({allele_name})' = {latex(simplified_next_freq)}",
            "conclusion": f"f({allele_name})' = {latex(original_freq)} = f({allele_name})"
        }
        
        next_gen_calcs.append(calculation_detail)
    
    derivation_steps.append({
        "step": 4,
        "title": "Next Generation Allele Frequency Recovery",
        "content": "Calculate next generation allele frequencies from current genotype frequencies:",
        "calculations": next_gen_calcs
    })
    
    # Step 5: Final conclusion
    derivation_steps.append({
        "step": 5,
        "title": "Conservation Conclusion",
        "content": "Since all allele frequencies remain unchanged:",
        "conclusion": [
            f"f({name})' = f({name}) = {latex(freq)}" for name, freq in zip(allele_names, alleles)
        ],
        "final": "Therefore, genotype frequencies in the next generation will be identical to the current generation, demonstrating Hardy-Weinberg equilibrium."
    })
    
    return {
        "derivation_steps": derivation_steps,
        "genotype_frequencies": {name: latex(freq) for name, freq in genotype_freqs.items()},
        "num_alleles": num_alleles
    }

# Example usage and testing
if __name__ == "__main__":
    # Test the function
    result = derive_hw_equilibrium_detailed(3)
    
    # Print formatted output
    for step in result["derivation_steps"]:
        print(f"\n{'='*50}")
        print(f"STEP {step['step']}: {step['title']}")
        print('='*50)
        print(step['content'])
        
        if 'constraint' in step:
            print(f"\n{step['constraint']}")
        
        if 'equations' in step:
            for eq in step['equations']:
                print(f"  {eq}")
        
        if 'calculation' in step:
            print(f"\n{step['calculation']}")
        
        if 'calculations' in step:
            for calc in step['calculations']:
                print(f"\nFor allele {calc['allele']}:")
                print(f"  {calc['formula']}")
                print(f"  Substituting:")
                for sub in calc['substitution']:
                    print(f"    {sub}")
                print(f"  {calc['expanded']}")
                print(f"  {calc['simplified']}")
                print(f"  ∴ {calc['conclusion']}")
        
        if 'conclusion' in step:
            for conc in step['conclusion']:
                print(f"  {conc}")
        
        if 'final' in step:
            print(f"\n{step['final']}")

@app.route("/derive", methods=["GET"])
def derive():
    num_alleles = int(request.args.get("num_alleles", 3))
    detailed = request.args.get("detailed", "true").lower() == "true"
    
    if detailed:
        # Use the new comprehensive derivation
        result = derive_hw_equilibrium_detailed(num_alleles)
        return jsonify({"derivation": result})
    else:
        # Use the old simple derivation
        result = derive_hw_equilibrium(num_alleles)
        return jsonify({"derivation": result})

def allele_check(freqs, num_alleles=3):
    """
    Recover allele freqs from genotype freqs.
    Works for any num_alleles.
    """
    next_freqs = []
    steps = {}
    
    for i in range(num_alleles):
        allele = f"A{i+1}"
        total = 0
        step_parts = []
        
        for g, f in freqs.items():
            if g == allele+allele:  # homozygote
                total += f
                step_parts.append(f"{g}={f:.6f}")
            elif allele in g:       # heterozygote
                total += 0.5 * f
                step_parts.append(f"½×{g}={0.5*f:.6f}")
        
        next_freqs.append(total)
        steps[f"{allele}_next_steps"] = " + ".join(step_parts)
    
    result = {f"A{i+1}'": val for i, val in enumerate(next_freqs)}
    result["steps"] = steps
    
    return result

# -------------------------------
# API Endpoints
# -------------------------------

@app.route("/api/genotypes")
def api_genotypes():
    """
    Example call:
    /api/genotypes?p=0.2&q=0.3&r=0.5
    """
    try:
        p = float(request.args.get("p", 0))
        q = float(request.args.get("q", 0))
        r = float(request.args.get("r", 0))
    except Exception:
        return jsonify({"error": "Invalid parameters"}), 400

    try:
        freqs = genotype_frequencies(p, q, r)
        
        # Create verbose genotype output
        verbose_genotypes = {}
        verbose_genotypes["A1A1"] = f"p² = {p:.3f}² = {freqs['A1A1']:.6f}"
        verbose_genotypes["A2A2"] = f"q² = {q:.3f}² = {freqs['A2A2']:.6f}"
        verbose_genotypes["A3A3"] = f"r² = {r:.3f}² = {freqs['A3A3']:.6f}"
        verbose_genotypes["A1A2"] = f"2pq = 2 * {p:.3f} * {q:.3f} = {freqs['A1A2']:.6f}"
        verbose_genotypes["A1A3"] = f"2pr = 2 * {p:.3f} * {r:.3f} = {freqs['A1A3']:.6f}"
        verbose_genotypes["A2A3"] = f"2qr = 2 * {q:.3f} * {r:.3f} = {freqs['A2A3']:.6f}"

        # Create verbose conservation output
        p_next = freqs['A1A1'] + 0.5*freqs['A1A2'] + 0.5*freqs['A1A3']
        q_next = freqs['A2A2'] + 0.5*freqs['A1A2'] + 0.5*freqs['A2A3']
        r_next = freqs['A3A3'] + 0.5*freqs['A1A3'] + 0.5*freqs['A2A3']
        
        conservation_steps = {}
        conservation_steps["p_next"] = {
            "formula": f"p' = p² + ½(2pq) + ½(2pr)",
            "substitution": f"p' = {freqs['A1A1']:.6f} + ½({freqs['A1A2']:.6f}) + ½({freqs['A1A3']:.6f})",
            "value": p_next
        }
        conservation_steps["q_next"] = {
            "formula": f"q' = q² + ½(2pq) + ½(2qr)",
            "substitution": f"q' = {freqs['A2A2']:.6f} + ½({freqs['A1A2']:.6f}) + ½({freqs['A2A3']:.6f})",
            "value": q_next
        }
        conservation_steps["r_next"] = {
            "formula": f"r' = r² + ½(2pr) + ½(2qr)",
            "substitution": f"r' = {freqs['A3A3']:.6f} + ½({freqs['A1A3']:.6f}) + ½({freqs['A2A3']:.6f})",
            "value": r_next
        }

    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    return jsonify({
        "input": {"p": p, "q": q, "r": r},
        "genotypes": freqs,
        "verbose_genotypes": verbose_genotypes,
        "conservation": conservation_steps
    })

@app.route("/api/grid")
def api_grid():
    genotype = request.args.get("genotype", "A1A1")
    n = int(request.args.get("n", 50))

    grid = genotype_grid(n=n, genotype=genotype)
    p_vals = np.array(grid["p_vals"]).tolist()
    q_vals = np.array(grid["q_vals"]).tolist()

    # Replace NaN/Inf with None before JSON export
    Z = np.array(grid["z_matrix"])
    Z = np.where(np.isfinite(Z), Z, None).tolist()

    # Flatten points for ternary plot
    points = []
    for i, p in enumerate(grid["p_vals"]):
        for j, q in enumerate(grid["q_vals"]):
            r = 1 - p - q
            freq = grid["z_matrix"][j][i]
            if r >= 0 and np.isfinite(freq):
                points.append({
                    "p": float(p),
                    "q": float(q),
                    "r": float(r),
                    "freq": float(freq)
                })

    return jsonify({
        "p_vals": p_vals,
        "q_vals": q_vals,
        "z_matrix": Z,
        "points": points
    })

@app.route("/api/grid_csv")
def api_grid_csv():
    genotype = request.args.get("genotype", "A1A1")
    n = int(request.args.get("n", 50))

    if genotype not in ["A1A1","A2A2","A3A3","A1A2","A1A3","A2A3"]:
        return jsonify({"error": "Invalid genotype"}), 400

    grid = genotype_grid(n=n, genotype=genotype)

    # Create CSV in memory
    def generate():
        header = ["p","q","r","freq"]
        yield ",".join(header) + "\n"
        for i, p in enumerate(grid["p_vals"]):
            for j, q in enumerate(grid["q_vals"]):
                r = 1 - p - q
                freq = grid["z_matrix"][j][i]
                if r >= 0 and np.isfinite(freq):
                    yield f"{p},{q},{r},{freq}\n"

    return Response(generate(), mimetype="text/csv")

if __name__ == "__main__":
    app.run(debug=True, port=2025)
