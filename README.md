# Partition-Binomial Divisibility Analyzer

This repository contains the Python code to computationally investigate the divisibility of the central binomial coefficient $\binom{2n}{n}$ by the integer partition function $p(n)$. The project is based on the paper:

### Pre-print (Zenodo) : [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17162713.svg)](https://doi.org/10.5281/zenodo.17162713)
* **DOI** - 10.5281/zenodo.17162713
* **URL** - https://doi.org/10.5281/zenodo.17162713

## Abstract

This project explores the number theory problem: For which non-negative integers $n$ does the partition function $p(n)$ divide the central binomial coefficient $\binom{2n}{n}$? The code implements a novel and efficient algorithm using partial factorization up to $2n$ and two rigorous obstruction tests (prime bound and q-adic valuation) to make a large-scale search computationally feasible. This script was used to verify that the only solutions for $n \le 20,000$ are the known twelve values in the set $\mathcal{S} = \{0, 1, 2, 4, 5, 6, 8, 11, 14, 17, 18, 26\}$.

## Features

-   **Efficient Divisibility Checking:** Implements a fast algorithm to test the condition $p(n) \mid \binom{2n}{n}$.
-   **Large-Scale Search:** Capable of analyzing the problem for a large range of $n$ (tested up to $n=20,000$).
-   **Resumable Computation:** Generates and stores partition data in chunked `.pkl` files, allowing computation to be stopped and resumed.
-   **Congruence Analysis:** Includes functions to compute the exact value of $\binom{2n}{n} \pmod{p(n)}$ for cases where all prime factors of $p(n)$ are known.
-   **Pattern Discovery:** Can be used to search for other patterns and congruences in the remainders.

## Repository Contents

-   `analyze_partitions.py`: The core Python script to run the computation and analysis.
-   `partition_data/`: Directory where the computed `.pkl` data chunks are stored (created on first run).
-   `README.md`: This documentation file.
-   `requirements.txt`: A list of the Python dependencies required to run the script.
-   `computation_log.txt`: The log file capturing the time taken to run the script.
-   `LICENSE`: The MIT License file for the project.

## Requirements

-   Python 3.8+
-   Libraries listed in `requirements.txt`: `pandas`, `sympy`, `numba`

## Installation

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/your-username/partition-binomial-divisibility-analyzer.git](https://github.com/your-username/partition-binomial-divisibility-analyzer.git)
    cd partition-binomial-divisibility-analyzer
    ```

2.  **Install the required dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## Usage

Run the script from the command line. It will create the `partition_data` directory, compute and save the data, and then print the final analysis.

```bash
python analyze_partitions.py
```


## Generating Computational Logs (in Google Colab)

The `computation_log.txt` file in this repository was generated to provide a verifiable record of the script's output. To reproduce this log, follow these steps in a Google Colab notebook:

1.  **Prepare the main cell.** In a single code cell, place the `%%capture` command at the very top, followed by the entire content of the `analyze_partitions.py` script.

    ```python
    %%capture output_log
    # The '%%capture' line MUST be the first line.

    # Paste the entire contents of analyze_partitions.py here...
    import time
    from pathlib import Path
    # ... etc.
    ```

2.  **Run the cell.** Execute this cell. It will run silently (producing no visible output) as it captures the entire process into the `output_log` variable. This may take a long time to complete.

3.  **Save the log.** In a new cell below the first one, run the following code to save the captured output to a file.

    ```python
    with open('computation_log.txt', 'w') as f:
        f.write(output_log.stdout)
    ```

This will create the `computation_log.txt` file in your Colab session, which you can then download.

## Citation

If you use this work, please cite the paper using the Zenodo archive.

@misc{naladiga_venkat_2025_17162713,
  author       = {Naladiga Venkat, Arvind},
  title        = {Divisibility of Central Binomial Coefficients by
                   the Partition Function: A Computational Study and
                   Evidence for Finiteness
                  },
  month        = sep,
  year         = 2025,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.17162713},
  url          = {https://doi.org/10.5281/zenodo.17162713},
}

---

## License

The content of this repository is dual-licensed:

- **MIT License** for `analyze_partitions.py` See the [LICENSE](LICENSE) file for details.
- **CC BY 4.0** (Creative Commons Attribution 4.0 International) for all other content (results.txt, README, etc.)



## Author

- **Arvind N. Venkat** - [arvind.venkat01@gmail.com](mailto:arvind.venkat01@gmail.com)
