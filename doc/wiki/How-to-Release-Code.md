## Major/Minor release
In this example, we release version `gamer-2.3.0`, with the previous version being `gamer-2.2.1`.

1. Create a new branch `v2.3.x`
2. Bump `VERSION` from `gamer-2.3.dev` to `gamer-2.4.dev` in `include/Macro.h` at the `main` branch
3. Bump `VERSION` from `gamer-2.3.dev` to `gamer-2.3.0` in `include/Macro.h` at the `v2.3.x` branch
4. Go to `Release`
[[/images/Release.png | alt=Release]]
5. Click `Draft a new release`
[[/images/DraftANewRelease.png | alt=DraftANewRelease]]
6. Create a new tag `gamer-2.3.0` and select `v2.3.x` as the target branch
[[/images/ChooseNewTag.png | alt=ChooseNewTag]]
7. Click `Generate release notes`; polish the content if needed
[[/images/GenerateReleaseNote_Major_Minor.png | alt=GenerateReleaseNote_Major_Minor]]
8. Add the release title `GAMER 2.3.0`
9. Click `Publish release`
[[/images/PublishRelease.png | alt=PublishRelease]]


## Maintenance release
In this example, we release version `gamer-2.2.2`, with the previous version being `gamer-2.2.1`.

1. Cherry-pick the necessary commits from the latest `main` branch into `v2.2.x`
2. Bump `VERSION` from `gamer-2.2.1` to `gamer-2.2.2` in `include/Macro.h` at the `v2.2.x` branch
3. Go to `Release`
[[/images/Release.png | alt=Release]]
4. Click `Draft a new release`
[[/images/DraftANewRelease.png | alt=DraftANewRelease]]
5. Create a new tag `gamer-2.2.2` and select `v2.2.x` as the target branch
[[/images/ChooseNewTag.png | alt=ChooseNewTag]]
6. Click `Generate release notes`; polish the content if needed
[[/images/GenerateReleaseNote_Maintenance.png | alt=GenerateReleaseNote_Maintenance]]
7. Add the release title `GAMER 2.2.2`
8. Click `Publish release`
[[/images/PublishRelease.png | alt=PublishRelease]]
