# README: Turbulent Boundary Layer Trailing-Edge Noise Calculation

## 1. General Information
* **Author:** SichengZ
* **Date:** December 28, 2025
* **Entry Point:** The main entry point for the code is `Pred_Green.m`.

---

## 2. Project Description
This code is designed to calculate the turbulent boundary layer trailing-edge (TE) noise. The calculation supports both straight and serrated trailing edges.

### **Serration Geometry**
The serration is described by the function:
$$x=hF(y/\lambda)$$
Where:
* $h$ is the half root-to-tip amplitude.
* $\lambda$ is the serration wavelength.
* $F(\eta)$ is the periodic shape function of the serration with a period of 1.

> **Note:** This code only works for **sawtooth serrations**. All notations used are consistent with the symbols used in the associated research paper.

---

## 3. Program Outputs
The program outputs include:
* The scattered noise Power Spectral Density (PSD), obtained by summing the spanwise modes and multiplying by the input turbulent boundary layer wavenumber-frequency spectrum of the flat plate.

---

## 4. Models and Sub-functions
The code implements the following models and functions:
* **Input Spectrum:** Chase's model is used for the input TBL spectrum, as defined in the function `Pi_Chase.m`.
* **Sub-functions:**
    * `g_n` corresponds to $g_n^{nj}$ in the paper.
    * `g_01` corresponds to $g_0^{(01)}$.
    * `g_02` corresponds to $g_0^{(02)}$.
    * `g_03` corresponds to $g_0^{(03)}$.
    * `g_04` corresponds to $g_0^{(04)}$.
    * `g_r`  represents the Curle integral result for the reflected wave.
    * `g_in` represents the Curle integral result for the incident wave.

---

## 5. Examples
* The example code `example.m` provides guidance for users.
