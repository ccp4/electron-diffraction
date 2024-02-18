import setuptools

with open("README.md", "r") as fh:
   long_description = fh.read()

setuptools.setup(
    name="ccp4ED",
    version="1.3",
    author="Tarik Ronan Drevon",
    author_email="tarik.drevon@stfc.ac.uk",
    description="Electron diffraction utilities",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://pypi.org/project/ccp4ED",
    project_urls={
        'Documentation': 'https://debloch.readthedocs.io/en/latest/',
        'Source':'https://github.com/ccp4/electron-diffraction',
    },

    # packages=['multislice','scattering','wallpp'],
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={'multislice/data':['splines.npy']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License ",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['wheel','numpy','scipy','matplotlib','colorama','pandas',
    'crystals','TarikDrevonUtils','easygui','tifffile','pickle5','bs4',
    'cbf','mrcfile','gemmi','opencv-python',
    ],
)
