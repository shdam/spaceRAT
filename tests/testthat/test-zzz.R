test_that(".onUnload unloads spaceRATScaffolds namespace", {
    # Load the package to ensure it's loaded
    library(spaceRAT)
    library(spaceRATScaffolds)

    # Call the .onUnload function
    .onUnload()

    # Check if the namespace is unloaded
    expect_error(detach("package:spaceRATScaffolds", unload = TRUE), "invalid 'name' argument")
})
