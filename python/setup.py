from setuptools import setup, find_packages

setup(
    name="universalfft",
    version="1.0.0",
    description="Boteler (2012) compliant FFT/IFFT wrappers for geoscience",
    author="UniversalFFT contributors",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.23",
    ],
    extras_require={
        "dev": [
            "pytest>=7",
            "pytest-cov>=4",
            "scipy>=1.9",
        ],
    },
)
