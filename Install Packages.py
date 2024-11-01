import subprocess
import sys

# Get the current Python version
python_version = sys.version
print(python_version)

def check_and_install_package(package_name):
    if package_name == 'sklearn':
        package_name_for_pip = 'scikit-learn'
    else:
        package_name_for_pip = package_name

    try:
        __import__(package_name)
        print(f"'{package_name}' is already installed.")
    except ModuleNotFoundError:
        print(f"'{package_name}' is not installed. Installing now...")
        try:
            if package_name == 'balanced_clustering':
                # Use --ignore-requires-python only for balanced_clustering
                subprocess.check_call([sys.executable, "-m", "pip", "install", package_name_for_pip, "--ignore-requires-python"])
            else:
                subprocess.check_call([sys.executable, "-m", "pip", "install", package_name_for_pip])

            print(f"'{package_name}' has been installed successfully.")
            __import__(package_name)
        except subprocess.CalledProcessError as e:
            print(f"Failed to install '{package_name}'. Error: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")

# Example usage: Check and install numpy
check_and_install_package('numpy')
check_and_install_package('sklearn')
check_and_install_package('balanced_clustering')