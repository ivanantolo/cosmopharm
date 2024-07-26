# Minimal Clone with Sparse Checkout, Filtered Blobs, and Shallow Submodules
Follow these steps to clone only the necessary parts of the repository and its submodules:

## Steps

1. **Clone the repository with minimal data:**

    ```
    git clone --filter=blob:none --depth 1 --shallow-submodules --no-checkout https://github.com/usnistgov/COSMOSAC.git
    ```

2. **Change to the repository directory:**

    ```
    cd COSMOSAC
    ```

3. **Enable sparse checkout:**

    ```
    git config core.sparseCheckout true
    ```

4. **Create and populate the `.git/info/sparse-checkout` file:**
```
src/
include/
externals/
tests/
setup.py
CMakeLists.txt
```

5. **Checkout the specified files:**

    ```
    git checkout
    ```

6. **Initialize and update submodules with shallow clone and filtered blobs:**

    ```
    git submodule update --init --recursive --depth 1 --filter=blob:none
    ```

By following these steps, you will clone only the specified files and directories from the repository and its submodules, minimizing the data transfer and storage requirements.
