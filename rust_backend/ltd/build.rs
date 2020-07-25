fn main() {
    println!("cargo:rustc-link-search=/home/hirschva/MG5/pynloop_MadNkLO/PLUGIN/pyNLoop/libraries/fjcore");
    println!("cargo:rustc-link-search=/home/hirschva/MG5/HEP_softs/scs/out");
    println!("cargo:rustc-link-search=/home/hirschva/MG5/HEP_softs/ecos");
    println!("cargo:rustc-link-search=/usr/lib/python3.7/config-3.7m-x86_64-linux-gnu");
    println!("cargo:rustc-link-lib=stdc++");
    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=gcc");
    println!("cargo:rustc-link-lib=gfortran");
    println!(r"cargo:rustc-link-search=/home/hirschva/tools/Cuba-4.2/lib");
}
