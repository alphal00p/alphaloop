fn main() {
    println!("cargo:rustc-link-search=../libraries/fjcore");
    println!("cargo:rustc-link-search=../libraries/ecos");
    println!("cargo:rustc-link-search=../libraries/scs/out");
    println!("cargo:rustc-link-search=../libraries/Cuba-4.2");
    println!("cargo:rustc-link-lib=stdc++");
    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=gcc");
    println!("cargo:rustc-link-lib=gfortran");
}
