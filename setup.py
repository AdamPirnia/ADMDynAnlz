# setup.py
from setuptools import setup
from Cython.Build import cythonize

setup(
  name="md_pipeline_ext",
  ext_modules=cythonize([
    "main_functions/coordinates_extract.pyx",
    "main_functions/unwrap_coords.pyx",
    "main_functions/COM_calc.pyx",
    "main_functions/alpha2_MSD.pyx",
  ], compiler_directives={"language_level": "3"}),
)

