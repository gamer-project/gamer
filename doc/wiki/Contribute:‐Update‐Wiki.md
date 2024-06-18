We are glad you have something new to contribute to GAMER, and you would like to teach people how to use it!

Before getting started, I recommend you have some basic knowledge of [repository, branch, fork](https://docs.github.com/en/repositories/creating-and-managing-repositories/about-repositories), [action](https://docs.github.com/en/actions), and [workflow](https://docs.github.com/en/actions/using-workflows).

The following is the outline of this document:
* [Setup](#Setup)
* [Edit Wiki](#Edit-Wiki)
* [Contribute](#Contribute)
* [Reference](#Reference)

## Setup
1. **Check GAMER version**
   * Please make sure your current GAMER commit is later than commit `13409ab33b12d84780076b6a9beb07317ca145f1` (committed on 2024/06/10).

1. **Create first page**
   * Click the Wiki page and create the first page

   [[images/CreateWiki.png]]

1. **Create a token for the action bot**
   * Go to `Setting` of your account > `Developer setting` > `Personal access tokens` > `Generate new token`
   * Please check the `repo` and the `workflow` options.
   * You might want to set the `Expiration` to `No expiration`.

   [[images/CreateToken.png]]

   * Remember to save the token since it will only be shown once!

1. **Create a repository secret token**
   * Go to `Setting` of the repository > `Security` > `Secrets and variables` > `Actions` > `Repository secrets`
   
   [[images/CreateSecret.png]]

   * Click `New repository secrets`, and then you will see the following. Please replace `<your personal access token>` to your token at `Secret` and make sure the `Name` of the token is `TOKEN_FOR_WIKI`. 

   [[images/SetToken.png]]

1. **Enable actions**
   * Enable workflows to run in the forked repository.
   * Click `Actions` > click the green button.
   ![image](https://github.com/ChunYen-Chen/CheckNode/assets/70311975/9e58d4a8-3248-4ceb-81ff-276a6943149d)

1. **Initialize wiki**
   * Click `Action` > `Workflows` > `Copy doc to wiki` > `Run workflow` > choose `Branch: master` > Click green `Run workflow`. Once the workflow is done, the wiki is also updated to the master branch version.

   [[images/InitializeWiki.png]]

## Edit Wiki
We provide two examples of editing the wiki pages: from GitHub website and through `gollum`. In the following examples, we call the new branch you would like to contribute as `new_contribution_branch`.
1. **From GitHub website**
   We treat the wiki page on GitHub as a local in this method, and we use the workflow to approach the `pull` and `push` behavior of the wiki.
   - **Copy wiki doc to wiki page**

     This step is like the `git checkout new_contribution_branch` but the GitHub wiki website version. 
     * Click `Action` > `Workflows` > `Copy doc to wiki` > `Run workflow` > choose `new_contribution_branch` branch > Click green `Run workflow`. Once the workflow is done, the wiki is also updated to the specific branch version.

      [[images/CopyDocFromNewBranch.png]]

   - **Edit on GitHub wiki website**
   - **Copy wiki to wiki doc**

     This step is like the `git push` to your branch.
     * Click `Action` > `Workflows` > `Sync wiki to doc` > `Run workflow` > choose `new_contribution_branch` branch > Click green `Run workflow`. Once the workflow is done, the docs in the specific branch are updated to the latest wiki.

      [[images/CopyWikiToNewBranch.png]]

1. **Through `gollum`**
   - Install [gollum](https://github.com/gollum/gollum).
   - Clone wiki git. You may find the wiki url at the bottom right of the Wiki page.
   - Edit by gollum.
     * `gollum /path/to/wiki`
     * Open `http://localhost:4567` in your browser.
   - Push your changes to Wiki git
   - **Copy wiki to wiki doc**

     This step is like the `git push` to your branch.
     * Click `Action` > `Workflows` > `Copy wiki to doc` > `Run workflow` > choose the specific branch > Click green `Run workflow`. Once the workflow is done, the docs in the specific branch are updated to the latest wiki.

     [[images/CopyWikiToNewBranch.png]]

## Contribute
   1. File a PR to the GAMER public version.

## Reference
1. [Bi-directional Wiki Sync Action](https://github.com/marketplace/actions/bi-directional-wiki-sync-action)
1. [Create Pull Requests For Your GitHub Wiki](https://nimblehq.co/blog/create-github-wiki-pull-request)
1. [gollum](https://github.com/gollum/gollum)