# 🌊 myRIOMAR_dev

**myRIOMAR_dev** is a Python toolkit designed to process satellite data to support marine research as part of the RiOMar project ([RiOMar website](https://riomar.lsce.ipsl.fr/)). 
It enables automated data download, validation, regional mapping, plume detection, time series decomposition.

<p align="center">
  <img src="https://img.shields.io/github/license/louis-terrats/myRIOMAR_dev" alt="License">
  <img src="https://img.shields.io/github/last-commit/louis-terrats/myRIOMAR_dev" alt="Last Commit">
</p>

---

## 🚀 Features

- 📥 **Data Downloading** – Automate download of satellite data  
- ✅ **Data Validation** – Quality control of satellite data using matchups with *insitu* measurements  
- 🗺️ **Regional Mapping** – Create and customize regional satellite maps  
- 🌫️ **Plume Detection** – Identify river plumes and their extent from remote sensing data  
- 📈 **X11 Analysis** – Apply time series decomposition (trend, seasonal, residual) using the X11 method ([reference here](https://www.sciencedirect.com/science/article/abs/pii/S0967063711000380))
- 📊 **Figure Generation** – Ready-to-publish figures for scientific reports and articles  

---

## 📦 Installation

### Using Conda (recommended)

#### For Unix/macOS or generic setup
```
conda env create -f myRIOMAR.yml
conda activate myRIOMAR
```

#### Using pip
Clone the repository and install dependencies:
```
git clone https://github.com/louis-terrats/myRIOMAR_dev.git
cd myRIOMAR_dev
pip install -r requirements.txt
```

---

## 🛠 Requirements

Make sure the following are installed:

- Python ≥ 3.8
- Packages listed in requirements.txt
- R (used for some post-processing steps, mainly plotting)

---

## 🧭 Project Structure

```
myRIOMAR_dev/
├── _0_data_downloading/       # Scripts to download satellite data
├── _1_data_validation/        # Quality control of satellite data using matchups with *insitu* measurements  
├── _2_regional_maps/          # Generate regional maps
├── _3_plume_detection/        # Detect river plumes and analyze morphology
├── _4_X11_analysis/           # Time series decomposition (X11 method)
├── _5_Figures_for_article/    # Reproducible figure generation for a publication led by David Doxaran.
├── _99_common/                # R and Python utilities common to all modules
├── main.py                    # Execution script of package functions
```

---

## 📚 Usage

All the work can be done with the 'main.py' script.

This script will be cleaned up in the coming weeks, to make it easier to use.

---

## 📄 License

This project is licensed under the MIT License.

---

## 👤 Author

Developed and maintained by Louis Terrats.

For scientific or technical questions, feel free to open an issue or reach out via GitHub.
