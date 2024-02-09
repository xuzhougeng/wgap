from os import path
from os.path import join
import setuptools

# read the README
this_directory = path.abspath(path.dirname(__file__))
with open(join( this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setuptools.setup(
    name="wgap", # Replace with your own username
    version="0.7.1",
    url="https://github.com/xuzhuogeng/wgap",
    project_urls={
        "Bug Tracker": "https://github.com/xuzhuogeng/wgap/issues"
    },
    author="Zhou-geng Xu",
    author_email="xuzhougeng@163.com",
    description="A whole genome annotation pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['wgap'],
    package_data = {'': ["wgap/*",]},
    data_files = [ (".", ["README.md"] ) ],
    include_package_data = True,
    python_requires=">=3.7",
    entry_points = {
        'console_scripts': [
            'wgap = wgap.wgap:cli'
        ]
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)
