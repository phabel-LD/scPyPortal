from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="scpyportal",
    version="0.1.0",
    author="Phabel Antonio Lopez Delgado",
    author_email="phabel2001@gmail.com | phabel@lcg.unam.mx",
    description="Interactive Single-Cell Analysis Platform built with Shiny for Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.10",
    install_requires=requirements,
    include_package_data=True,
    package_data={
        "": ["static/*", "static/**/*", "*.html", "*.css", "*.js"],
    },
    entry_points={
        "console_scripts": [
            "scpyportal=app:main",
        ],
    },
    keywords="single-cell, bioinformatics, shiny, python, scanpy, interactive",
    url="https://github.com/yourusername/scpyportal",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/scpyportal/issues",
        "Source": "https://github.com/yourusername/scpyportal",
    },
)