include("../src/FakeModule.jl")
include("../src/FakeModule2.jl")

using Main.FakeModule
using Main.FakeModule2

FakeModule.f(3)

f(2)

g(2)

h(2)
