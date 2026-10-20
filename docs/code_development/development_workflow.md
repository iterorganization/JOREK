---
title: "Development Workflow"
nav_order: 1
parent: "Code Development"
layout: default
render_with_liquid: false
---


# Main Development Workflow

## Committing Guidelines
- See also [**Coding Guidelines**](develop#coding-guidelines)

<!--- - [Slides by Guido (2016-04)](assets/develop//huijsmans_jorek_general_meeting_2016.pdf)
- [Slides by Matthias (2020-04)](assets/develop/2020-04-mhoelzl-development-coordination.pdf) with some issues coming up regularly and suggested solutions --->


## Repository structure
- The [main repository](https://github.com/iterorganization/JOREK) is public and serves as the central source for JOREK.
- Developers will create forks of the main repository for their individual development work.
- Collaborators can be added to forks to enable team-based development on specific features or bug fixes.
- **Developments in the forks need to be merged back into the main repository as soon as possible**.

---

## Setting up the workflow
1. **Create a github account**
    - If you do not have one already, setup a github account and make sure your username contains your last name and initial.

2. **Fork the Main Repository**:
   - If you don't have a fork already, create a new one under your user account. 
   - Navigate to the [main repository on GitHub](https://github.com/iterorganization/JOREK.git).
   - Click the "Fork" button to create a personal copy of the repository under your GitHub account.
   - Add collaborators to your fork if you are working closely with someone else.

3. **Clone Your Fork**:
   - Clone your fork to your local machine:

```bash
     git clone https://github.com/<your-username>/<repository-name>.git
     cd <repository-name>
```

4. **Set Up the documentation**:
   - During your development, you will contribute to the documentation. You can preview it as it will appear in 
   the main repository by setting 'Settings'->'Pages'->'Source' to 'GitHub Actions' *in your fork*. 
   For faster iteration, however, we recommend building the documentation locally.
   You can find instructions for this in the [documentation setup guide](/docs/contributing.md).

5. **Set Up the Main Repository as an Upstream Remote**:
   - Add the main repository as an upstream remote to keep your fork in sync:

```bash
     git remote add upstream https://github.com/iterorganization/JOREK.git
```

6. **Sync Your Fork with the Main Repository**:
   - Regularly fetch and merge changes from the main repository:

```bash
     git fetch upstream
     git merge upstream/develop
```


## Development and Bugfixes

1. **Develop Features or Fix Bugs**:
   - Create a new branch for each feature or bug fix:

```bash
     git checkout -b <feature|bugfix>/branch-name
```

   - Commit your changes to the branch:

```bash
     git add .
     git commit -m "Description of changes"
```

2. **Push Changes to Your Fork**:
   - Push your branch to your fork:

```bash
     git push origin <feature-or-bugfix-branch>
```

3. **Create a Pull Request**:
   - Open a pull request (PR) from your fork's branch to the main repository's `develop` branch.
       - From the main page of your fork's repository a button with "Contribute" will appear
       - Click on the button and then "Open Pull Request"
       - Fill out the Pull Request Form describing the Pull Request in a way that is understandable to the whole community
   - Resolve potential merge conflicts and verify the automatic regression tests are passing.
   - Ensure the PR is reviewed and approved before merging by at least two people. 


## Best practices
- Merge changes from the main repository into your fork frequently to avoid large merge conflicts. This is especially important for long developments
- Keep feature and bug fix branches small and focused.
- Ensure all changes are reviewed (by at leat 2 people) before merging into the main repository.
- Make sure to create new [regression tests](nrt) for your developments.
- Bug fixes should be developed independently of feature development.
- Merge bug fixes from the main repository into your fork as soon as they are available:

```bash
    git fetch upstream
    git merge upstream/develop
```

## Using the GitHub Issue Tracker
Use the [GitHub Issue Tracker](https://github.com/iterorganization/JOREK/issues) to bring up important bugfixes and discussions on new developments or to make a feature request.

1. **Creating Issues**:
  - Provide a clear and concise title and description.
  - Include relevant details such as steps to reproduce a bug, expected behavior, or the motivation for a new feature.
  - Mention specific users with `@username` to make them aware of the issue.
  - Use appropriate labels (e.g., `bug`, `enhancement`, `discussion`) to categorize issues.

2. **Linking Pull Requests**:
  - Reference issues in pull requests using keywords like `Fixes #<issue-number>` to automatically close the issue when the PR is merged.
  - Ensure the issue is resolved and all discussions are addressed before closing.
  - More infos here about [linking an issue to a PR](https://docs.github.com/en/issues/tracking-your-work-with-issues/using-issues/linking-a-pull-request-to-an-issue)




## Migration to GitHub
During the migration to GitHub, the history of commit hashes will change. Follow these steps to ensure a smooth transition:

1. **Push Your Branches Before Migration (end of July)**:
   - Before the migration, push all your local branches to the current repository:

```bash
     git push origin <branch-name>
```

2. **Retrieve Your Branch After Migration**:
   - After the migration, contact **jorek-devel@jorek.eu** to retrieve your new branch from the private GitHub repository where old branches are saved during the migration period, **please include your GitHub username in the email**.

3. **Set Up Your Fork with the Old Branch**:
   - Be sure to have followed the steps in [Setting Up the Workflow](#setting-up-the-workflow)

   - Clone your fork of the new main repository:

```bash
     git clone git@github.com:<your-username>/<fork-repo>.git
     cd <fork-repo>
```

   - Add the private repository (containing the old branches) as a remote:

```bash
     git remote add old-repo git@github.com:N-Schwarz/jorek-backup.git
```

   - Fetch your old branch:

```bash
     git fetch old-repo <branch-name>
```

   - Create a branch in your fork based on the old branch:

```bash
     git checkout -b <branch-name> old-repo/<branch-name>
```

   - Push the branch to your fork:

```bash
     git push -u origin <branch-name>
```

  - You might need to use the following command to make the branch track your fork

```bash             
     git branch feature/m600refluid  -u origin/feature/m600refluid
```

   The old branch and the public repository share the same history, so the
   commits attach where they always were: nothing is rewritten and there are
   no conflicts to resolve. 

---


