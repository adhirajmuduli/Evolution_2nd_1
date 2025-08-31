Hardy-Weinberg Equilibrium Explorer
Overview

The Hardy-Weinberg Equilibrium Explorer is an interactive web application designed to demonstrate and visualize the core principles of population genetics. This tool allows users to explore the Hardy–Weinberg equilibrium for both two- and three-allele systems, providing a clear, step-by-step breakdown of genotype frequency calculations, conservation of allele frequencies, and dynamic data visualization.

This project is an ideal educational resource for students and researchers in biology, genetics, and bioinformatics. It demystifies complex genetic concepts by offering a hands-on experience with:

    Derivation Display: See a complete, LaTeX-formatted algebraic derivation of the Hardy–Weinberg principle for multi-allele systems.

    Numeric Verification: Input custom allele frequencies (p, q, r) and see the calculated genotype frequencies and a verification that allele frequencies are conserved across generations.

    Data Visualization: Plot interactive 3D surfaces and 2D ternary plots to visualize the relationship between allele frequencies and genotype frequencies.

Features

    Interactive Controls: Easily switch between two- and three-allele systems and adjust allele frequencies with intuitive input fields.

    Live Calculations: Instantly calculate genotype frequencies and verify the conservation of allele frequencies as you adjust the inputs. The detailed calculations are displayed in a human-readable format, making the underlying math transparent.

    Dynamic Visualizations: The application generates interactive plots that visually represent the Hardy–Weinberg principle. The 3D surface plot shows how a genotype's frequency changes with varying allele frequencies, while the 2D ternary plot visualizes the same data on a population simplex.

    Responsive Design: The interface is optimized for both desktop and mobile devices, ensuring a seamless experience for all users.

    Dark/Light Mode: A user-controlled theme toggle allows you to switch between a light and dark interface to suit your viewing preference.

Technology Stack

The Hardy-Weinberg Equilibrium Explorer is built on a robust and modern technology stack to ensure a fast, reliable, and interactive experience.

    Frontend:

        HTML5: For structuring the web page.

        CSS3: For styling and creating a clean, professional user interface. Custom CSS with variables enables the dark/light mode functionality.

        JavaScript: Powers all the interactive elements, including input validation, live calculations, and data fetching.

        Plotly.js: A powerful, open-source JavaScript graphing library used to create the interactive 3D and 2D visualizations.

        MathJax: An industry-standard JavaScript display engine for rendering LaTeX-formatted mathematical equations directly in the browser.

    Backend:

        Flask: A lightweight Python web framework that handles API requests, processes calculations, and serves the frontend files.

        Numpy: Used for efficient numerical calculations and to generate the data for the 3D plots.

Getting Started

To run this project locally, follow these steps:

    Clone the repository:
    git clone https://github.com/your-username/hardy-weinberg-explorer.git

    Navigate to the project directory:
    cd hardy-weinberg-explorer

    Install the required Python packages:
    pip install -r requirements.txt

    Run the Flask application:
    python ass1.py

The application will be accessible at http://127.0.0.1:2025 in your web browser.
