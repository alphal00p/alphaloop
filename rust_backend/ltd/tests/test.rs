// Use ltd
extern crate ltd;

#[cfg(test)]
mod test_arbitrary_module_name {
    // import relevant functions and variables 
    use ltd::topologies::LoopLine;

    #[test]
    fn test_name() {
        // very complicated test
        assert_eq!(2 + 2, 4);
    }
}
