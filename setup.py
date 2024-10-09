from setuptools import setup

setup(name='MassChargeCuration',
      version='0.9',
      description='Library to be used for mass and charge balancing a model',
      #url='',
      #author='',
      #author_email='',
      license='MIT',
      packages=['MCC.DataCollection', 'MCC.DataCollection.Requests', 'MCC.Balancing', 'MCC.ModelInterface', 'MCC.ReportGeneration', 'MCC'],
      install_requires=['numpy', 'pandas', 'requests', 'matplotlib', 'z3-solver', 'python-libsbml', 'tqdm'],
      zip_safe=False)
