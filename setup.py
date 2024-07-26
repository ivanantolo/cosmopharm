from setuptools import setup, find_packages

# Dynamically set the long_description from readme_pypi.md, if available
try:
    with open('./docs/README_pypi.md', encoding='utf-8') as f:
        long_description = f.read()
except FileNotFoundError:
    # Fallback to a default description or leave as is if setup.cfg has enough info
    long_description = None

if __name__ == "__main__":
    setup(
        long_description=long_description,
        long_description_content_type='text/markdown',
        packages=find_packages(where="src"),
        package_dir={"": "src"},
    )
