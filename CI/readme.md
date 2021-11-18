# CI - Scripts

These scripts provide all CI functions.
They also allow developers to run them localy or the provided dockerfile

**All scripts have to be run from the repository path!**

## Environment

To be able to check all functionalities, a virtual python environment `env-CI` will be created.

```sh
./CI/run-env.sh
```

## Checking Code Format

The `C++` formatation is checked via clang-format.
The style is define in the `.clang-format` file n the repository.

For the python formation `flake8` is used.

```sh
./CI/run-code-check.sh
```

## Checking Build

The `fastpli` package will be installed inside the `env-CI` environment.

```sh
./CI/run-build.sh
```

## Checking Tests

All tests inside the `.\tests\` path will be run.

```sh
./CI/run-pytest.sh
```

## Checking Examples

All examples inside the `./examples` path will be checked for runnability.

```sh
./CI/run-examples.sh
```

## API and WIKI:

- make docs
- make wiki
- switch to ssh credential
- push new version
- notebook md tutorials have to be formated by `pretty`
