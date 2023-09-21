# General considerations
- When compiling into a single file, accessory files are extracted to a temp location before running. This location changes each run, and the app doesn't know beforehand where. So to use accessory files you need to include in main.py:
  - import os, sys
  - from kivy.resources import resource_add_path, resource_find
  - Under the `if __name__ == '__main__'` block, include the lines
    - `if hasattr(sys, '_MEIPASS'):
        resource_add_path(os.path.join(sys._MEIPASS))`
- When working with data files, we need to include them into the compiled file. Since the script to do this is general and includes everything recursively, it makes it easier to have a `Code` and a `Build` folder in the top-most level. `Code` is the venv, and `Build` is where all of the compilation happens. Otherwise the inclusion code copies a second version of all python libraries (desired or not) into the exe.
- The program icon was drawn in icon.svg, saved as icon.png, then converted to icon.ico using https://www.convertico.com/

# To compile on Windows
- Set up virtual environment `python -m venv Code`, cd in, `Scripts\activate.bat`.
- Install dependencies `pip install kivy`, `pip install --upgrade pyinstaller`.
- Transfer relevant files.
  - Under Code: main.py, icon.ico.
  - Create Code\app_data, transfer: aln_key_layout.kv.
  - Create Code\app_scripts, transfer: align.py, BLOSUM.py, sequ.py
  - Create Code\dev_files, transfer: icon.ico

- Create `Build` folder under the root, and cd in.
- Generate the spec file: `python -m PyInstaller --onefile --name AlnKey --windowed --icon C:\Users\curra\compilation\aln_key\Code\dev_files\icon.ico C:\Users\curra\compilation\aln_key\Code\main.py`
- Modify the file `AlnKey.spec`:
  - At the top add `from kivy_deps import sdl2, glew`
  - Within the Analysis() block, add to the hiddenimports list:
    - 'win32timezone' needed for some reason to allow filechooser to work
  - Within the EXE block, add the 2 lines:
  exe = EXE(
      pyz,
      `Tree('C:\\Users\\curra\\compilation\\aln_key\\Code\\app_data\\', 'app_data'),`
      ...,
      [],
      `*[Tree(p) for p in (sdl2.dep_bins + glew.dep_bins)],`
      ...,
  )
  - Without the 2nd 'app_data' arg in the Tree() call above, all files would be siblings to main.py in the compiled hierarchy (once it is unzipped at runtime). This preserves keeping those files in a child folder called 'app_data'.

- To reduce the final exe file size (only saves ~1.5MB):
  - At the top add `from kivy.tools.packaging.pyinstaller_hooks import get_deps_minimal, get_deps_all, hookspath, runtime_hooks`
  - In the `Analysis()` block, delete the lines: `binaries=`, `hiddenimports=`, `excludes=`.
  - At the end of the kwargs in the `Analysis()` block, add `**get_deps_minimal(video=None, audio=None, camera=None, spelling=None)`, which provides the 3 entries just deleted.

- Compile everything with: `python -m PyInstaller --clean AlnKey.spec`.
  - Without the `clean` option, the build would hang forever on `INFO: Building PKG (CArchive) AlnKey.pkg` and the subdirectory `build` would keep growing forever.
- The exportable exe is found in the `dist` folder.
