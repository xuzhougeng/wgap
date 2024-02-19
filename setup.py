from os import path
from src.wgap import __version__
from setuptools import find_packages, setup


# Read the README for the long description
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="wgap",
    version= __version__,
    url="https://github.com/xuzhuogeng/wgap",
    project_urls={
        "Bug Tracker": "https://github.com/xuzhuogeng/wgap/issues"
    },
    author="Zhou-geng Xu",
    author_email="xuzhougeng@163.com",
    description="A whole genome annotation pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_dir={"": "src"},  # Specify the root directory for packages
    packages=find_packages(where="src"),  # Find packages in the src directory
    package_data={
        "": ["*.md", "envs/*", "rules/*.smk", "*.yaml", "*.csv"],  # Include additional files
        "wgap": ["*.py"],  # Include all Python files in the wgap package
    },
    include_package_data=True,
    python_requires=">=3.7",
    entry_points={
        'console_scripts': [
            'wgap = main:cli'
        ]
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)
