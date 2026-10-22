use std::sync::Arc;

use crate::document::{Block, Category, Document};

/// An owned reference to one block of a shared [`Document`].
///
/// The handle holds the document through an [`Arc`], so the block stays
/// readable for as long as the handle exists.
///
/// # Examples
///
/// ```
/// use std::sync::Arc;
///
/// use mmcif_cif::{BlockHandle, Document};
///
/// let document = Arc::new(Document::parse("data_a\n_x.y 1\ndata_b\n_x.y 2\n")?);
/// let blocks = BlockHandle::all(&document);
/// drop(document);
///
/// assert_eq!(blocks.len(), 2);
/// assert_eq!(blocks[1].view().name(), "b");
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Debug, Clone)]
pub struct BlockHandle {
    document: Arc<Document>,
    index: usize,
}

/// An owned reference to one category of a block in a shared [`Document`].
///
/// # Examples
///
/// ```
/// use std::sync::Arc;
///
/// use mmcif_cif::{BlockHandle, Document};
///
/// let document = Arc::new(Document::parse("data_a\nloop_\n_atom.id\n1\n2\n")?);
/// let block = BlockHandle::new(document, 0).unwrap();
/// let atoms = block.category("ATOM").unwrap();
///
/// assert_eq!(atoms.name(), "atom");
/// assert_eq!(atoms.view().row_count(), 2);
/// # Ok::<(), mmcif_cif::ParseError>(())
/// ```
#[derive(Debug, Clone)]
pub struct CategoryHandle {
    block: BlockHandle,
    query: Box<str>,
}

impl BlockHandle {
    /// Returns a handle to the block at `index` in `document`, or `None` if
    /// `index` is out of range.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::sync::Arc;
    ///
    /// use mmcif_cif::{BlockHandle, Document};
    ///
    /// let document = Arc::new(Document::parse("data_a\n_x.y 1\n")?);
    ///
    /// assert_eq!(BlockHandle::new(Arc::clone(&document), 0).unwrap().index(), 0);
    /// assert!(BlockHandle::new(document, 1).is_none());
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn new(document: Arc<Document>, index: usize) -> Option<Self> {
        document.block(index)?;
        Some(BlockHandle { document, index })
    }

    /// Returns a handle to each block of `document`, in order.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::sync::Arc;
    ///
    /// use mmcif_cif::{BlockHandle, Document};
    ///
    /// let document = Arc::new(Document::parse("data_a\n_x.y 1\ndata_b\n_x.y 2\n")?);
    /// let names: Vec<_> = BlockHandle::all(&document)
    ///     .iter()
    ///     .map(|block| block.view().name().to_owned())
    ///     .collect();
    ///
    /// assert_eq!(names, ["a", "b"]);
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn all(document: &Arc<Document>) -> Vec<BlockHandle> {
        (0..document.len())
            .map(|index| BlockHandle {
                document: Arc::clone(document),
                index,
            })
            .collect()
    }

    /// Returns the document that holds the block.
    pub fn document(&self) -> &Arc<Document> {
        &self.document
    }

    /// Returns the position of the block in its document.
    pub fn index(&self) -> usize {
        self.index
    }

    /// Returns a view of the block.
    pub fn view(&self) -> Block<'_> {
        self.document
            .block(self.index)
            .expect("block handles are created for existing blocks")
    }

    /// Returns a handle to the category called `name`, or `None` if the block
    /// has no such category.
    ///
    /// `name` is matched as by [`Block::category`].
    ///
    /// # Examples
    ///
    /// ```
    /// use std::sync::Arc;
    ///
    /// use mmcif_cif::{BlockHandle, Document};
    ///
    /// let document = Arc::new(Document::parse("data_a\n_entry.id 1ABC\n")?);
    /// let block = BlockHandle::new(document, 0).unwrap();
    ///
    /// assert_eq!(block.category("_entry").unwrap().name(), "entry");
    /// assert!(block.category("cell").is_none());
    /// # Ok::<(), mmcif_cif::ParseError>(())
    /// ```
    pub fn category(&self, name: &str) -> Option<CategoryHandle> {
        self.view().category(name)?;
        Some(CategoryHandle {
            block: self.clone(),
            query: name.into(),
        })
    }
}

impl CategoryHandle {
    /// Returns the handle of the block that holds the category.
    pub fn block(&self) -> &BlockHandle {
        &self.block
    }

    /// Returns the category name, as [`Category::name`] gives it.
    pub fn name(&self) -> &str {
        self.view().name()
    }

    /// Returns a view of the category.
    pub fn view(&self) -> Category<'_> {
        self.block
            .view()
            .category(&self.query)
            .expect("documents are immutable, so a query that resolved once resolves again")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn handles_keep_the_document_alive() {
        let document =
            Arc::new(Document::parse("data_a\n_x.y 1\ndata_b\nloop_\n_z.w\n1\n2\n").unwrap());
        let blocks = BlockHandle::all(&document);
        drop(document);
        assert_eq!(blocks.len(), 2);
        let category = blocks[1].category("Z").unwrap();
        assert_eq!(category.name(), "z");
        assert_eq!(category.view().row_count(), 2);
        assert_eq!(blocks[0].view().name(), "a");
        assert!(blocks[0].category("z").is_none());
        assert!(BlockHandle::new(Arc::clone(blocks[0].document()), 2).is_none());
    }

    #[test]
    fn categories_whose_names_start_with_an_underscore_resolve() {
        let document = Arc::new(
            Document::parse(
                "data_a
__x.y 1
_x.z 2
",
            )
            .unwrap(),
        );
        let block = BlockHandle::new(document, 0).unwrap();
        assert_eq!(block.view().category_names(), ["_x", "x"]);
        let odd = block.category("_x").unwrap();
        assert_eq!((odd.name(), odd.view().items()), ("_x", vec!["y"]));
        let plain = block.category("x").unwrap();
        assert_eq!((plain.name(), plain.view().items()), ("x", vec!["z"]));
        assert_eq!(
            block.category("_x_missing").map(|c| c.name().to_owned()),
            None
        );
    }
}
