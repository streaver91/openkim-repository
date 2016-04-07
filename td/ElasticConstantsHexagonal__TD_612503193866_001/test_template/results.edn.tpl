[
{
    "property-id" "tag:staff@noreply.openkim.org,2014-04-15:property/bulk-modulus-isothermal-hexagonal-crystal-npt"
    "instance-id" 1
    
    "short-name" {
        "source-value"  ["@<crystal_structure>@"]
    }
    "species" {
        "source-value"  ["@<species>@"]
    }
    "a" {
        "source-value"  @<lattice_constant_a>@
        "source-unit"   "angstrom"
    }
    "c" {
        "source-value"  @<lattice_constant_c>@
        "source-unit"   "angstrom"
    }
    "basis-atom-coordinates" {
        "source-value"  @<basis_coordinates>@
    }
    "space-group"  {
        "source-value"  "@<space_group>@"
    }
    "temperature" {
        "source-value"  0
        "source-unit"  "K"
    }
    "cauchy-stress"  {
        "source-value"  [0 0 0 0 0 0]
        "source-unit"   "GPa"
    }
    "isothermal-bulk-modulus" {
        "source-value"  @<B>@
        "source-unit"   "@<units>@"
        "source-std-uncert-value"  @<B_sig>@
    }
}
]
