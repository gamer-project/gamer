We are glad you have something new to contribute to GAMER, and you would like to teach people how to use it!

Before getting started, I recommend you have some basic knowledge of [repository, branch, fork](https://docs.github.com/en/repositories/creating-and-managing-repositories/about-repositories), [action](https://docs.github.com/en/actions), and [workflow](https://docs.github.com/en/actions/using-workflows).

The followings are the outline of this document:
* [Introduction](#Introduction)
* [Setup](#Setup)
* [Edit Wiki](#Edit-Wiki)
* [Contribute](#Contribute)
* [Reference](#Reference)

## Introduction
In the following document, the Wiki stands for the Wiki page on the GitHub website, and the doc stands for the markdown files under `docs/wiki`. Remember the doc and the Wiki are NOT always synchronized, you will have to use the action to copy either from doc to Wiki or from Wiki to doc. 

## Setup
**_This setup only needs to be done once (the first time) !!!_**

1. **Check GAMER version**
   * Please make sure you have the latest GAMER.

1. **Create the first page**
   * Click the Wiki page and create the first page

   [[images/CreateWiki.png]]

1. **Create a token for the action bot**
   * Go to `Setting` of your account > `Developer setting` > `Personal access tokens` > `Generate new token`
      - NOTE: We use the `Tokens(classic)` in this example.
   * Please check the `repo` and the `workflow` options.
      - NOTE: You can set the `Note` of this token freely
   * You might want to set the `Expiration` to `No expiration`.
   * Click the green `Generate token` at the bottom.

   [[images/CreateToken.png]]

   * Remember to save the token since it will only be shown once!

1. **Create the repository secret token and email**
   * Go to `Setting` of your forked gamer repository > `Security` > `Secrets and variables` > `Actions` > `Repository secrets`
   
   [[images/CreateSecret.png]]

   * Click `New repository secrets`, and then you will see the following. Please replace `<your_personal_access_token>` with the token generated in the previous step at `Secret` and make sure the `Name` of the secret is `TOKEN_FOR_WIKI`. 

   [[images/SetToken.png]]

   * Click `New repository secrets`. Please replace `<your_email_address>` with your account email address at `Secret` and make sure the `Name` of the secret is `MY_EMAIL`.
      - NOTE: We will not know your address since it only exists in your repository. This step is only for recording your contribution by the action.

   [[images/SetMail.png]]

1. **Enable actions (workflows)**
   * Click `Actions` > click the green button.
   ![image](https://github.com/ChunYen-Chen/CheckNode/assets/70311975/9e58d4a8-3248-4ceb-81ff-276a6943149d)

1. **Initialize wiki**
   * Click `Action` > `Workflows` > `Copy doc to wiki` > `Run workflow` > choose `Branch: main` > Click green `Run workflow`. Once the workflow is done, the wiki is also updated to the `main` branch version.

   [[images/InitializeWiki.png]]

## Edit Wiki
We provide three examples of editing the wiki pages: through `gollum` (recommended), from the GitHub website, and on the local terminal. In the following examples, we call the new branch you would like to contribute as `new_contribution_branch`.

1. **Through `gollum` (recommended)**
   - Install [gollum](https://github.com/gollum/gollum).
   - Click `Action` > `Workflows` > `Copy doc to wiki` > `Run workflow` > choose `new_contribution_branch` branch > Click green `Run workflow`. Once the workflow is done, Wiki is also updated to the `new_contribution_branch` branch version.
   - Clone your forked Wiki git. You may find the Wiki URL at the bottom right of the Wiki page.

     [[images/WikiGitLocation.png]]

   - Edit by `gollum`. 
     * `gollum /path/to/wiki`
     * Open `http://localhost:4567` in your browser.
   - Push your changes to the forked Wiki git
   - **Copy wiki to doc**

     This step is like `git push` to your branch.
     * Click `Action` > `Workflows` > `Copy wiki to doc` > `Run workflow` > choose `new_contribution_branch` branch > Click green `Run workflow`. Once the workflow is done, the docs in the specific branch are updated to the latest Wiki.

     [[images/CopyWikiToNewBranch.png]]

1. **From GitHub website**

   We treat the Wiki page on GitHub as a local in this method, and we use the workflow to approach the `pull` and `push` behavior of the Wiki.
   - **Copy wiki doc to wiki page**

     This step is like `git checkout new_contribution_branch` but the GitHub Wiki website version. 
     * Click `Action` > `Workflows` > `Copy doc to wiki` > `Run workflow` > choose `new_contribution_branch` branch > Click green `Run workflow`. Once the workflow is done, the Wiki is also updated to the `new_contribution_branch` branch version.

      [[images/CopyDocFromNewBranch.png]]

   - **Edit on GitHub Wiki website**

      See [About wikis - GitHub Docs](https://docs.github.com/en/communities/documenting-your-project-with-wikis/about-wikis) for more details.
   - **Copy wiki to doc**

     This step is like `git push` to your branch.
     * Click `Action` > `Workflows` > `Sync wiki to doc` > `Run workflow` > choose `new_contribution_branch` branch > Click green `Run workflow`. Once the workflow is done, the docs in the `new_contribution_branch` branch are updated to the latest wiki.

      [[images/CopyWikiToNewBranch.png]]

1. **On the local terminal**

   Just like how you edit the code on the local terminal. The files of Wiki can be found at `docs/wiki`. Once you push your changes to `origin`, you may preview your changes on the Wiki page by doing:
   * Click `Action` > `Workflows` > `Copy doc to wiki` > `Run workflow` > choose `new_contribution_branch` branch > Click green `Run workflow`. Once the workflow is done, the Wiki is also updated to the `new_contribution_branch` branch version.

## Contribute
   Before you file a new PR, please check that all the hyperlinks work, the images are clear enough, and all the related pages about the PR have been updated.
   1. File a PR to the GAMER public version.
      * NOTE: The contribution is made through the files in the doc (e.g. `doc/wiki/`), not through the Wiki page of the forked repository.

## Reference
1. [Bi-directional Wiki Sync Action](https://github.com/marketplace/actions/bi-directional-wiki-sync-action)
1. [Create Pull Requests For Your GitHub Wiki](https://nimblehq.co/blog/create-github-wiki-pull-request)
1. [gollum](https://github.com/gollum/gollum)