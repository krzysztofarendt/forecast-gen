import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "fcastgen",
    version = "0.0.1dev",
    author = "Krzysztof Arendt",
    author_email = "krza@mmmi.sdu.dk",
    description = "Forecast generation functions",
    long_description = long_description,
    url = "https://github.com/krzysztofarendt/forecast-gen",
    packages = setuptools.find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
)