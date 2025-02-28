fn main() {
    if cfg!(feature = "python_api") {
        pyo3_build_config::add_extension_module_link_args();
    }
    println!("cargo:rustc-link-search=native=/opt/local/Library/Frameworks/Python.framework/Versions/3.11/lib");
    println!("cargo:rustc-link-lib=python3.11");
    println!("cargo:rustc-link-search=../libraries/fjcore");
    println!("cargo:rustc-link-search=../libraries/ecos");
    println!("cargo:rustc-link-search=../libraries/scs/out");
    println!("cargo:rustc-link-search=../libraries/Cuba-4.2.1");
    println!("cargo:rustc-link-lib=stdc++");
    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=gcc");
    println!("cargo:rustc-link-lib=gfortran");
    println!("cargo:rustc-link-lib=gcc_s");
}
