# AlnKey
Software to visualize sequence alignments and perform many common tasks

- This program was written using Kivy v2.1.0

## Installation
This repo is not designed to be cloned for general use. To install, download a pre-compiled binary from the distributions at XX

## Developer setup
- Make new root folder, probably named `aln_key`, cd in.
- Clone the github files into a new subfolder named Code `git clone https://github.com/dave-the-scientist/AlnKey.git Code`
- Set it up as virtual environment `python -m venv Code`.
- cd in and activate venv `Scripts\activate.bat` (Windows)
- Install dependencies `pip install kivy`, `pip install --upgrade pyinstaller`.
	- There used to be a bug in pyinstaller preventing execution of kivy programs at least up until v5.13, which required the command to be `pip install "pyinstaller==5.6.2"`. Bug appears to be fixed as of at least v6.8.