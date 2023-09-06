from setuptools import setup, find_packages

setup(
    name="CageCavityCalc",
    version="1.0.2",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "scikit-learn",
        "rdkit"
    ],
    entry_points={
        "console_scripts": [
            "CageCavityCalc = CageCavityCalc.__main__:main"
        ]
    }
)