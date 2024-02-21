# Network Coherence Project

**Authors:**
[Alex Fuster](https://github.com/Alex-Fuster)
[Katherine Hébert](https://github.com/katherinehebert)
[Francis Banville](https://github.com/FrancisBanville)
[Dominique Caron](https://github.com/DominiqueCaron)
[Dominique Gravel](https://github.com/DominiqueGravel)
[Gracielle Higino](https://github.com/graciellehigino)
[Pablo Silva](https://github.com/tutudocerrado)

The objective of our study is to develop a metric that quantifies the degree of association between a community response to a specific environment and its species interactions. The primary goal we envision for this metric is to serve as an indicator of the predictability and stability of communities in the face of environmental changes. Additionally, the metric could also allow for exploring numerous relevant questions on species' responses to the biotic and abiotic environments.

You can read more about our theoretical approach [here](supporting_material/network_coherence_12feb24.Rmd)

## Folder Structure

The project follows the following folder structure:


- `a_code/`: Contains the source code files for our analyses.
 - `functions/`: Contains the stand-alone functions we call in our scripts.
- `b_data/`: Contains the data files for the project. This folder is not watched by git. Do not upload data into our repository. 
    - `raw/`: Only raw data - do not manipulate! 
    - `clean/`: Only cleaned data (outputs of raw data manipulations).
- `c_outputs/`: Contains all the outputs of our analyses, including figures.
- `d_manuscript/`: Contains the manuscript files (source and rendered), including *.bib file and templates.
- `e_supporting_material/`: Contains any extra supporting material for reference.

Additional files can be found in our [Google Drive](https://drive.google.com/drive/folders/1pmDwl0QXEWIqchEwgF7LnLsoww7ru1QV?usp=sharing), such as raw data, reference papers, the program of our in-person meeting, and technical guides.


## Naming Convention

The project follows the following naming convention:

- Always use lowercase

- Code files: Use descriptive names that reflect the purpose of the file. Enumerate according to the order in which the scripts should be run. Ex.: "00-download_data.R", "01-clean_occurrence_data.R", etc.

- Folders: only first-level folders get an alphabetical order prefix (e.g., "a_ code/"). Subfolders should be ordered as needed.


**1. Case Sensitivity:**
   - Use lowercase letters for all file and folder names to ensure consistency and avoid potential issues with case sensitivity across different platforms and file systems.
   - Use hyphens (-) to separate categories of information on your file name, and underscores (_) to separate words in the same category (e.g., in the "01-clean_occurrence_data.R" file name, 01 is a prefix intended to ordinate the code, while "clean_occurrence_data" is a description of what the code does. 

**2. Folder Structure:**
   - Organize folders logically to reflect the hierarchy and structure of the software.
   - Use descriptive folder names that clearly indicate the contents or purpose of each folder.
   - Avoid nesting folders too deeply to maintain clarity and ease of navigation. 

**3. Naming Folders:**
   - Begin folder names with a noun or a descriptive term that indicates the contents.
   - Keep folder names concise yet descriptive to convey their purpose at a glance.

**4. Naming Files Differently:**
   - Use unique and descriptive names for files to distinguish them from others.
   - Incorporate relevant keywords or terms related to the file's content or functionality.
   - Include version numbers or dates in file names if necessary to differentiate between multiple versions or revisions.

**5. Separating Code Files:**
   - Clearly distinguish between different types of code files based on their purpose or programming language.
   - Use consistent file extensions to denote the programming language (e.g., .py for Python, .R for R, .js for JavaScript).
   - Group related code files together within designated folders to facilitate easier maintenance and collaboration.

**Example File Naming:**

- `species_occurrence_data.csv`
- `cleaned_occ_data.csv`
- `02-correlation_analysis.R`
- `user_guide.pdf`
- `formating_reference.md`

Adhering to this naming convention will ensure consistency, clarity, and ease of use across the development team, facilitating efficient collaboration and maintenance of our project.


## Contributing

To contribute to this project, follow the rules delineated in the [CONTRIBUTING.md](CONTRIBUTING.md) file.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
