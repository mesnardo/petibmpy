# How to contribute to petibmpy

Welcome to the developer's guide of petibmpy!

## Adding new features and fixing bugs

All new features and bug fixes must go through a pull-request review procedure.
If you want to contribute to petibmpy, please fork the main [petibmpy](https://github.com/mesnardo/petibmpy) repository, make your changes on your fork, and then open a pull-request.

For new features and minor bugs (with small impact), the base branch of the pull-request should be the `develop` branch of the main repository.
(The `develop` branch will be merged into the `master` one once we are ready for a new release of petibmpy.)

For major bugs, the base branch should be the `master` branch of the main repository; it will be considered as a hotfix (bugfix) and a new version of petibmpy will be released as soon as possible by the maintainers with the micro number incremented.

New features should come with some kind of test or example to verify and/or validate the implementation.

## Reporting bugs and requesting new features

To report bugs, request new features, or simply ask questions, please open a GitHub issue on the main repository.

## Writing documentation

New classes, methods, and functions must be documented with doctrings.

You should also add code documentation whenever necessary; it will greatly help other developers to review your new features and bug fixes.
