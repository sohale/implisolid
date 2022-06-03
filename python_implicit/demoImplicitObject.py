from implicitObject import Heart

# Init new obj
new_obj = Heart()

# Marching Cubes it
new_obj.build_marching_cubes()

# Show it!

new_obj.show(representation="surface")
# or
new_obj.show() #which uses wireframe representation as a default.
