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

## 📁 Project Structure

myRIOMAR_dev/                                             
├── CODES/                           
│   ├── _0_data_downloading/         
│   ├── _1_data_validation/          
│   ├── _2_regional_maps/            
│   ├── _3_plume_detection/          
│   ├── _4_X11_analysis/             
│   └── _5_Figures_for_article/     
├── CONFIG/                          
│   ├── general.yaml                 
│   ├── paths.yaml                   
│   └── multithreading.yaml          
├── .Rhistory                        
├── .gitattributes                   
├── .gitignore                       
├── Environment.yml                  
├── LICENSE                          
├── README.md                        
├── __main__.py                      
├── requirements.txt                 
└── setup.py                         

---


## 📦 Installation

### 1) Clone the repository

```
git clone https://github.com/louis-terrats/myRIOMAR_dev.git
cd myRIOMAR_dev
```

### 2) Install dependencies

#### Using Conda (recommended)
```
conda env create -f Environment.yml
conda activate myRIOMAR
```

#### Using pip
```
pip install -r requirements.txt
```

---

## 🛠 Requirements

Make sure the following are installed:

- Python ≥ 3.8
- Packages listed in requirements.txt
- R (used for some post-processing steps, mainly plotting)

---

## 📚 Usage

### ⚙️ Configuration
All parameters live in CONFIG/*.yaml:

- general.yaml: data sources, sensors, variable lists, date range

- paths.yaml: Package_root, Data_save_dir, Results_dir

- multithreading.yaml: nb_of_cores_to_use, plume thresholds

_Edit these files to customize your workflow._

### 🎯 Command-Line Interface
After installation, run:

```
python -m myRIOMAR_dev <command> [options]
```
| Command      | Description                                       | Options       |
| ------------ | ------------------------------------------------- | ------------- |
| `download`   | Download satellite data                           | `--overwrite` |
| `validate`   | Match satellite with in-situ measurements         | `--redo`      |
| `maps`       | Generate regional maps                            | `--overwrite` |
| `plume`      | Detect plumes (dynamic or fixed threshold)        | `--dynamic`   |
| `timeseries` | Plot plume area time-series                       | *(none)*      |
| `x11`        | Apply X11 seasonal decomposition                  | *(none)*      |
| `figures`    | Generate figures for publication (Fig. 1…Fig. 10) | *(none)*      |

_Examples_:

```
python -m myRIOMAR_dev download --overwrite
python -m myRIOMAR_dev maps --overwrite
python -m myRIOMAR_dev plume --dynamic
```

---

## 📄 License

This project is licensed under the MIT License.

---

## 👤 Author

Developed and maintained by Louis Terrats.

For scientific or technical questions, feel free to open an issue or reach out via GitHub.
