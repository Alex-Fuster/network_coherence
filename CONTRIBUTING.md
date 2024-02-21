**Contributing Guidelines**

*Welcome to the Contributing Guidelines for the Network Coherence Working Group!*

This document outlines the best practices and guidelines for contributing code, documentation, and other enhancements to our project.

**Branching Strategy:**

1. **Main Branch (`main`):**
   - The `main` branch represents the latest stable version of the software.
   - Direct commits to the `main` branch are restricted. All changes must be made through pull requests.

2. **Feature Branches:**
   - Create a new branch for each new feature or enhancement you are working on.
   - Name feature branches descriptively, using lowercase letters and hyphens or underscores to separate words according to our naming convention (e.g., `documentation-new_analysis_method`).

3. **Bugfix Branches:**
   - For bug fixes, create a new branch from `main` named appropriately (e.g., `bugfix-issue-123`).
   - Ensure the branch name includes a reference to the related issue or bug number.

**Pull Request Guidelines:**

1. **Naming Convention:**
   - Name your pull requests clearly and concisely, summarizing the changes made (e.g., "Add spatial analysis functionality").

2. **Description:**
   - Provide a detailed description of the changes introduced by the pull request.
   - Include any relevant context, such as why the changes were made and how they address any associated issues.

3. **Reviewers:**
   - Assign at least one reviewer to your pull request. It can be your teammate!
   - Choose reviewers who are knowledgeable about the affected areas of the codebase.


4. **Documentation:**
   - Update any affected documentation to reflect the changes made by your pull request.
   - This includes code comments, README files, and user guides.

5. **Merging:**
   - Once approved by the reviewers and all checks pass, the pull request can be merged into `main`.
   - Squash commits if necessary to maintain a clean and organized commit history.

**General Guidelines:**

1. **Code Style:**
   - Follow the existing code style and conventions used in the project.
   - Use meaningful variable names and adhere to best practices for readability and maintainability ([see below](#Code-Style-Guidelines)).

2. **Commit Messages:**
   - Write clear and descriptive commit messages that explain the purpose of each change.
   - You can use [emojis](https://gitmoji.dev/) to help us quickly locate your changes by categories.
   - Start the commit message with a present-tense verb (e.g., "Add", "Fix", "Update").

3. **Respect and Professionalism:**
   - Treat all contributors with respect and professionalism.
   - Provide constructive feedback and be open to receiving feedback from others.


## Code Style Guidelines

Our research software project in R follows specific code style guidelines to ensure consistency, readability, and maintainability of the codebase. Adhering to these guidelines facilitates collaboration among developers and promotes a professional coding environment. Here are the key aspects of our R code style:

1. **Indentation:**
   - Use 2 spaces for each level of indentation. Avoid using tabs.

2. **Naming Conventions:**
   - Variable names should be descriptive and use underscores to separate words (e.g., `occurrence_data`, `correlation_index`).
   - Function names should be lowercase and use underscores to separate words (e.g., `process_data`, `generate_plots`).

3. **Comments:**
   - Include comments to explain complex logic, algorithms, or any non-obvious parts of the code. More importantly, explain ***WHY*** you made these choices.
   - Comments should be concise and written in complete sentences to enhance understanding.

4. **Whitespace:**
   - Use a single blank line to separate logical sections within functions and scripts.
   - Maintain consistent spacing around operators and after commas (e.g., `x <- 10`, not `x<-10`).

5. **Package Imports:**
   - Organize package imports at the top of the script or function.
   - Use `library()` or `require()` statements to load required packages.
   - Import functions from packages explicitly to improve code readability.

6. **Line Length:**
   - Limit lines to a maximum of 80 characters to ensure readability in most code editors.
   - Break long lines into multiple lines when necessary for clarity, using parentheses for line continuation.

7. **Consistency:**
   - Follow the existing code style used within the project to maintain uniformity throughout the codebase.
   - Adhere to any project-specific conventions or patterns established for R code.

By following these guidelines, we will avoid unnecessary merge conflicts and we will maintain a high level of quality and collaboration. Thank you for your contributions! :heart"\: