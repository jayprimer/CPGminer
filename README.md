![CPGminer: Complete Prokaryote Genome (CPG) Metadata Dashboard](/images/CPGlogo2.png)

# CPGminer: Complete Prokaryote Genome (CPG) Metadata Dashboard
CPGminer is a Python application that provides an interactive dashboard for Complete Prokaryote Genome (CPG) metadata. It's designed to offer an easy way to inspect and analyze genome features and taxonomy, promoting greater understanding and facilitating improved data-driven decisions in genome research and education.

## Getting Started
These instructions will get you a copy of the project up and running on your local machine.

### Prerequisites
- Python (version 3.7 or above)
- pip (Python package installer)

## Clone the repository
Clone the project to your machine:

```bash
git clone git@github.com:jayprimer/CPGminer.git
```

## Setting Up a Virtual Environment
To isolate your project and manage its dependencies, it's recommended to create a virtual environment.

### Create a virtual environment
Navigate to the project directory and run the following command to create a new virtual environment:

```bash 
python -m venv ./venv
```

### Activate the virtual environment
Activate the virtual environment using the following command:

For linux,
```bash
source ./venv/bin/activate
```

For windows,
```
.\venv\Scripts\activate
```

### Deactivate the virtual environment
When you're done working on the project, you can deactivate the virtual environment:
```bash
deactivate
```


## Installation
Once you have cloned the repository and activated your virtual environment, navigate to the project directory and install the required dependencies using pip:

```bash
pip install -r requirements.txt
```

## Running the Application
This application uses Streamlit to provide a user-friendly interface. After the installation of dependencies, you can run the application with the following command:

```bash
streamlit run main.py
```

This command starts the Streamlit server and you can interact with the application by opening your web browser to http://localhost:8501.

## License
This project is licensed under the MIT License - see the LICENSE.md file for details.
