GAMER is distributed on the
[GAMER GitHub repository](https://github.com/gamer-project/gamer).
To download the code, you will need to install
[git](https://help.github.com/articles/set-up-git)
(if not installed yet) and then type

```bash
git clone https://github.com/gamer-project/gamer
```

By default, you will be using the development version after `git clone`.
If you would like to switch to the stable version, type

```bash
git checkout stable
```

If you have downloaded the code previously and want to update
to the latest development or stable version, move to the code directory and type
one of the following commands.

```bash
# for the development version
git pull origin main

# for the stable version
git pull origin stable
```

See [guides.github.com](https://guides.github.com/) for guides to git and GitHub.

> [!CAUTION]
> To validate whether the downloaded code correctly retains the symbolic links
used by GAMER, try
> ```bash
> ls -l gamer/src/SelfGravity/GPU_Gravity/
> ```
> Correct: `CUPOT_ExternalAcc.cu -> ../CPU_Gravity/CPU_ExternalAcc.cpp` <br/>
> Incorrect: `CUPOT_ExternalAcc.cu`
>
> Try updating `git` if you get incorrect results. Adding either `--config core.symlinks=true`
or `--config core.symlinks=false` to your `git clone` command may also help.
