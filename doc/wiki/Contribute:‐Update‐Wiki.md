We are glad you have something new to contribute to GAMER, and you would like to teach people how to use it!

Before getting started, I recommend you have some basic knowledge of [repository, branch, fork](https://docs.github.com/en/repositories/creating-and-managing-repositories/about-repositories), [action](https://docs.github.com/en/actions), and [workflow](https://docs.github.com/en/actions/using-workflows).

The following is the outline of this document:
* [Setup]()
* [Edit Wiki]()
* [Contribute]()
* [Reference]()

## Setup
1. **Check GAMER version**
   * Please make sure your GAMER commit is later than commit `13409ab33b12d84780076b6a9beb07317ca145f1` (committed on 2024/06/10).

1. **Create first page**
   * Click the Wiki page and create the first page

   ![image](https://github.com/ChunYen-Chen/CheckNode/assets/70311975/c48243fb-d603-4475-9949-3a088561c5c1)

1. **Create a token for the action bot**
   * Go to `Setting` of your account > `Developer setting` > `Personal access tokens` > `Generate new token`
   * Please check the `repo` and the `workflow` options.
   * You might want to set the `Expiration` to `No expiration`.

   ![image](https://github.com/ChunYen-Chen/CheckNode/assets/70311975/5e19015d-5ebb-46b6-9a7c-eb3fff298527)

   * Remember to save the token since it will only be shown once!

1. **Create a repository secret token**
   * Go to `Setting` of the repository > `Security` > `Secrets and variables` > `Actions` > `Repository secrets`
   
   ![image](https://github.com/ChunYen-Chen/CheckNode/assets/70311975/1e7a7a4e-924f-442d-8f96-d48e4f1dc783)

   * Click `New repository secrets`, and then you will see the following. Please replace `<your personal access token>` to your token at `Secret` and make sure the `Name` of the token is `TOKEN_FOR_WIKI`. 

   <img width="1208" alt="Screenshot 2024-04-26 at 23 50 38" src="https://github.com/ChunYen-Chen/CheckNode/assets/70311975/adb7656b-ab65-40f5-8b25-37d302cf4e77">

1. **Enable actions**
   * Enable workflows to run in the forked repository.
   * Click `Actions` > click the green button.
   ![image](https://github.com/ChunYen-Chen/CheckNode/assets/70311975/9e58d4a8-3248-4ceb-81ff-276a6943149d)

1. **Initialize wiki**
   * Click `Action` > `Workflows` > `Copy doc to wiki` > `Run workflow` > choose `Branch: master` > Click green `Run workflow`. Once the workflow is done, the wiki is also updated to the master branch version.

   <img width="1414" alt="Screenshot 2024-04-27 at 11 20 57" src="https://github.com/ChunYen-Chen/CheckNode/assets/70311975/2cd3d25b-160c-4fd7-8d4f-359675ed99ee">

## Edit Wiki
We provide two examples of editing the wiki pages: from GitHub website and through `gollum`.
1. **From GitHub website**
   We treat the wiki page on GitHub as a local in this method, and we use the workflow to approach the `pull` and `push` behavior of the wiki.
   - **Copy wiki doc to wiki page**

     This step is like the `git checkout new_contrbution_branch` but the GitHub wiki website version. 
     * Click `Action` > `Workflows` > `Copy doc to wiki` > `Run workflow` > choose the specific branch > Click green `Run workflow`. Once the workflow is done, the wiki is also updated to the specific branch version.

      <img width="1414" alt="Screenshot 2024-04-27 at 11 20 57" src="https://github.com/ChunYen-Chen/CheckNode/assets/70311975/2cd3d25b-160c-4fd7-8d4f-359675ed99ee">

   - **Edit on GitHub wiki website**
   - **Copy wiki to wiki doc**

     This step is like the `git push` to your branch.
     * Click `Action` > `Workflows` > `Sync wiki to doc` > `Run workflow` > choose the specific branch > Click green `Run workflow`. Once the workflow is done, the docs in the specific branch are updated to the latest wiki.

      <img width="1416" alt="Screenshot 2024-04-27 at 11 22 00" src="https://github.com/ChunYen-Chen/CheckNode/assets/70311975/1d5482a2-8c33-41ee-bd2b-f7119924db81">

1. **Through `gollum`**
   - **Install [gollum](https://github.com/gollum/gollum)
   - clone wiki git
   - edit by gollum
   - push wiki git
   - **Copy wiki to wiki doc**

     This step is like the `git push` to your branch.
     * Click `Action` > `Workflows` > `Copy wiki to doc` > `Run workflow` > choose the specific branch > Click green `Run workflow`. Once the workflow is done, the docs in the specific branch are updated to the latest wiki.

      <img width="1416" alt="Screenshot 2024-04-27 at 11 22 00" src="https://github.com/ChunYen-Chen/CheckNode/assets/70311975/1d5482a2-8c33-41ee-bd2b-f7119924db81">

## Contribute
   1. file a PR

## Reference
1. [Bi-directional Wiki Sync Action](https://github.com/marketplace/actions/bi-directional-wiki-sync-action)
1. [Create Pull Requests For Your GitHub Wiki](https://nimblehq.co/blog/create-github-wiki-pull-request)
1. [gollum](https://github.com/gollum/gollum)