#!/usr/bin/env python

from pyroms_toolbox import ocean_in

def main():
#    fname = raw_input("Name of ocean.in file:")
    ocean_1 = ocean_in("ocean_palau1.in")
#    jname = raw_input("Name of json file:")
    ocean_1.write_json_dict("foo.json")
    ocean_1.write_json_dict("foo2.json", indent=2)

    ocean_1.merge_dicts([ocean_1, ocean_1])
    ocean_1.write_ocean_in('ocean_foo2.in')

def test_two():
    import subprocess
    ocean_1 = ocean_in('ocean_palau1.in')
    ocean_2 = ocean_in('ocean_palau2.in')

    ocean_1.merge_dicts([ocean_2])
    ocean_1.write_ocean_in('ocean_palau_1_2.in')

if __name__ == "__main__":
    test_two()
